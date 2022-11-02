library(magrittr)
library(dplyr)
library(Matrix)
ori_par <- par(no.readonly = TRUE)
# Set the application-level cache
shinyOptions(cache = cachem::cache_mem(max_size = 8 * 1024^3, evict = "lru"))

path_to_snakemake_folder <- ""
output_folder <- file.path(path_to_snakemake_folder, "results", "App")

# Define server ----
server <- function(input, output) {

    SS_file_path <-
        file.path(
            path_to_snakemake_folder,
            "results/Single_species_analysis/"
        )
    SS_ortho_path <-
        file.path(
            path_to_snakemake_folder,
            "results/Ortholog_analysis/"
        )
    MS_file_path <-
        file.path(
            path_to_snakemake_folder,
            "results/Multi_species_analysis/"
        )
    SS_files <- list.files(SS_file_path)
    MS_files <- list.files(MS_file_path)
    monosamples <- SS_files %>%
        grep(pattern = "cellranger_reanalysis_", value = TRUE) %>%
        sub(pattern = "cellranger_reanalysis_", replacement = "")
    multisamples <-
        list.files(file.path(MS_file_path, "all_plotting_tables")) %>%
        grep(pattern = "onlyUMAP", value = TRUE) %>%
        sub(pattern = "plotting_onlyUMAP_", replacement = "") %>%
        sub(pattern = ".csv", replacement = "")

    output$multisample_selection <- renderUI(
        if (input$UMAP_type == "multiple samples analysis (Seurat CCA)") {
            selectInput(
                inputId = "multisample",
                label = "Choose the multisample group:",
                choices = multisamples,
                selected = NULL
            )
        }
    )

    working_candidate_samples <- reactive({
        if (input$UMAP_type == "multiple samples analysis (Seurat CCA)") {
            df <- read.csv(
                file.path(
                    MS_file_path,
                    "all_plotting_tables",
                    paste0("plotting_", input$multisample, ".csv")
                )
            )
            out <- unique(df$Sample)
        } else {
            out <- monosamples
        }
        return(out)
    })

    output$plotting_sample_selection <- renderUI(
        if (
            input$UMAP_type == "multiple samples analysis (Seurat CCA)" &
            input$UMI_type == "ortholog UMI counts"
        ) {
            checkboxGroupInput(
                inputId = "plotting_sample",
                label = "Choose the sample(s):",
                choices = working_candidate_samples(),
                selected = working_candidate_samples()[1]
            )
        } else {
            radioButtons(
                inputId = "plotting_sample",
                label = "Choose the sample:",
                choices = working_candidate_samples(),
                selected = working_candidate_samples()[1]
            )
        }
    )

    working_feature_ids <- reactive({
        if (is.null(input$plotting_sample)) {
            out <- NULL
        } else {
            path_prefix <- ifelse(
                input$UMI_type == "ortholog UMI counts",
                file.path(SS_ortho_path, "all_UMI_tables", "orthologUMI_"),
                file.path(SS_file_path, "all_UMI_tables", "geneUMI_")
            )
            UMI_csv <- paste0(path_prefix, input$plotting_sample, ".csv")
            feature_id_list <-
                sapply(
                    UMI_csv,
                    function(x) {
                        con <- file(x, open = "r")
                        line <- readLines(con, n = 1) %>%
                            sub(pattern = "Barcode,", replacement = "") %>%
                            strsplit(split = ",") %>%
                            magrittr::extract2(1)
                        close(con)
                        return(line)
                    }
                )
            out <- unique(unlist(feature_id_list))
        }
        return(out)
    })

    rv <- reactiveValues(previous_plotting_feature = "")
    output$plotting_feature_selection <- renderUI(
        textInput(
            inputId = "plotting_feature",
            label = "Type a feature ID (gene ID or ortholog ID):",
            value = rv$previous_plotting_feature
        )
    )
    observe({
        rv$previous_plotting_feature <- input$plotting_feature
    }) %>%
        bindEvent(input$load_UMI_data)

    output$plotting_feature_selection_eg <- renderUI(
        helpText(paste("e.g. ", working_feature_ids()[1]))
    )

    # Get Sample_Barcode, is_plot, UMAP.1, UMAP.2
    get_plotting_df <- function(
        UMAP_type = input$UMAP_type,
        multisample = input$multisample,
        plotting_sample = input$plotting_sample
    ) {
        if (UMAP_type == "multiple samples analysis (Seurat CCA)") {
            plotting_csv <- file.path(
                MS_file_path,
                "all_plotting_tables",
                paste0("plotting_", multisample, ".csv")
            )
            plotting_df <- read.csv(plotting_csv)
            plotting_df$Sample_Barcode <- paste0(
                plotting_df$Sample, "_", plotting_df$Barcode
            )
            plotting_df$is_plot <-
                (plotting_df$Sample %in% plotting_sample)
        } else {
            plotting_csv <- file.path(
                SS_file_path,
                "all_plotting_tables",
                paste0("plotting_", plotting_sample, ".csv")
            )
            plotting_df <- read.csv(plotting_csv)
            plotting_df$Sample_Barcode <- paste0(
                plotting_sample, "_", plotting_df$Barcode
            )
            plotting_df$is_plot <- TRUE
        }
        return(plotting_df)
    }

    # Get Sample_Barcode, norm_factor
    get_norm_factor_df <- function(
        plotting_sample = input$plotting_sample
    ) {
        norm_factor_df_list <- sapply(
            simplify = FALSE,
            plotting_sample,
            function(x) {
                norm_factor_df <- read.csv(
                    file.path(
                        SS_file_path,
                        "all_norm_factor_tables", 
                        paste0("norm_factor_", x, ".csv")
                    )
                )
                norm_factor_df$Sample_Barcode <-
                    paste0(x, "_", norm_factor_df$Barcode)
                return(norm_factor_df)
            }
        )
        out <- bind_rows(norm_factor_df_list)
        return(out)
    }

    # Get UMI matrix (Sample_Barcode * feature_id)
    get_UMI_Matrix <- function(
        UMI_type = input$UMI_type,
        plotting_sample = input$plotting_sample
    ) {
        path_prefix <- ifelse(
            UMI_type == "ortholog UMI counts",
            file.path(SS_ortho_path, "all_UMI_tables", "orthologUMI_"),
            file.path(SS_file_path, "all_UMI_tables", "geneUMI_")
        )
        withProgress(
            message = "Loading UMI matrix:",
            value = 0,
            {
                UMI_df_list <- list()
                for (psi in seq_along(plotting_sample)) {
                    incProgress(
                        (2 * psi - 1) / (2 * length(plotting_sample)),
                        detail = plotting_sample[psi]
                    )

                    UMI_df <- read.csv(
                        paste0(path_prefix, plotting_sample[psi], ".csv")
                    )
                    UMI_df$Barcode <-
                        paste0(plotting_sample[psi], "_", UMI_df$Barcode)
                    UMI_df_list[[psi]] <- UMI_df

                    incProgress(
                        psi / length(plotting_sample),
                        detail = plotting_sample[psi]
                    )
                }
            }
        )
        all_UMI_df <- bind_rows(UMI_df_list)
        all_UMI_Matrix <- Matrix(
            as.matrix(all_UMI_df[, -1]),
            dimnames = list(
                all_UMI_df$Barcode,
                colnames(all_UMI_df)[-1]
            )
        )
        return(all_UMI_Matrix)
    }

    plot_figure_out <- function(
        plotting_df,
        norm_factor_df,
        UMI_Matrix,
        plotting_feature,
        color_range,
        main,
        output_clear,
        output_path
    ) {
        # Make color vector
        pp1 <- color_range[1]
        pp2 <- color_range[2]
        COL <- c(
            colorRampPalette(c("#C9D7EF", "#eeeeec"))(pp1),
            colorRampPalette(c("#eeeeec", "#ea4335"))(pp2 - pp1),
            colorRampPalette(c("#ea4335", "#ba1306"))(100 - pp2)
        )
        # Original blue: #4285f4

        # Reorder and extract data
        plotting_df <- plotting_df[plotting_df$is_plot, ]
        ref_Sample_Barcode <- plotting_df$Sample_Barcode
        stopifnot(ref_Sample_Barcode %in% norm_factor_df$Sample_Barcode)
        norm_factor_vector <- with(
            norm_factor_df,
            norm_factor[match(ref_Sample_Barcode, Sample_Barcode)]
        )
        if (any(is.na(norm_factor_vector))) {
            stop("Barcodes in norm_factor_csv do not match with plotting_csv.")
        }
        stopifnot(nrow(plotting_df) == length(norm_factor_vector))

        png(
            output_path,
            pointsize = 10, res = 600,
            width = 20, height = 15, units = "cm"
        )

        if (output_clear) {
            par(mai = c(0, 0, 0, 0))
        } else {
            par(mai = ori_par$mai)
        }

        if (plotting_feature %in% working_feature_ids()) {

            UMI_vector <- UMI_Matrix[ref_Sample_Barcode, plotting_feature]
            plotting_UMI_vector <- log2(UMI_vector * norm_factor_vector + 1)


            if (max(plotting_UMI_vector) == 0) {
                COL_ind <- rep(1, length(plotting_UMI_vector))
            } else {
                COL_ind <- round(
                    plotting_UMI_vector / max(plotting_UMI_vector) * 99 + 1,
                    digits = 0
                )
            }
            order_ind <- order(plotting_UMI_vector)
            plot(
                type = "n",
                plotting_df$UMAP.1,
                plotting_df$UMAP.2,
                xlab = "UMAP_1", ylab = "UMAP_2",
                axes = !output_clear
            )
            title(ifelse(output_clear, "", main), cex.main = 1.5, adj = 0)
            points(
                plotting_df$UMAP.1[order_ind],
                plotting_df$UMAP.2[order_ind],
                col = COL[COL_ind][order_ind],
                pch = 20, cex = 1
            )
            xlim <- par("usr")[1:2]
            ylim <- par("usr")[3:4]
            for (i in seq_along(COL)) {
                x <- seq(
                    qunif(0.65, xlim[1], xlim[2]),
                    qunif(0.95, xlim[1], xlim[2]),
                    length.out = length(COL) + 1
                )[c(i, i + 1)]
                y <- rep(ylim[2] + (ylim[2] - ylim[1]) * 7 / 100, 2)
                lines(x, y, lwd = 15, lend = "butt", col = COL[i], xpd = TRUE)
                if (i == pp1) lines(x, y, lwd = 15, lend = "butt", xpd = TRUE)
                if (i == pp2) lines(x, y, lwd = 15, lend = "butt", xpd = TRUE)
            }
            text(qunif(0.80, xlim[1], xlim[2]),
                ylim[2] + (ylim[2] - ylim[1]) * 2.5 / 100,
                "log2((UMI/totalUMI)*10000+1)",
                cex = 1, xpd = TRUE
            )
            text(qunif(0.95, xlim[1], xlim[2]),
                ylim[2] + (ylim[2] - ylim[1]) * 12 / 100,
                sprintf("%.2f", max(plotting_UMI_vector)),
                cex = 1.5, xpd = TRUE
            )
            text(qunif(0.65, xlim[1], xlim[2]),
                ylim[2] + (ylim[2] - ylim[1]) * 12 / 100,
                sprintf("%.2f", 0),
                cex = 1.5, xpd = TRUE
            )

        } else {

            plot(
                type = "n",
                plotting_df$UMAP.1,
                plotting_df$UMAP.2,
                xlab = "UMAP_1", ylab = "UMAP_2",
                axes = !output_clear
            )
            title(ifelse(output_clear, "", main), cex.main = 1.5, adj = 0)
            points(
                plotting_df$UMAP.1,
                plotting_df$UMAP.2,
                col = "#9f9f9f",
                pch = 20, cex = 1
            )
            legend("center", legend = "NA", cex = 4, bty = "n")

        }

        par(mai = ori_par$mai)
        dev.off()
    }

    working_UMI_Matrix <- reactive(get_UMI_Matrix()) %>%
        bindCache(input$UMI_type, input$plotting_sample) %>%
        bindEvent(input$load_UMI_data)

    output$UMI_data_in_use <- renderPrint({
        if (grepl("Matrix", class(working_UMI_Matrix()))) {
            out <-
                list(
                    `Type of UMAP` = input$UMAP_type,
                    `Type of UMI counts` = input$UMI_type,
                    `Plotting sample` = input$plotting_sample
                )
        }
        return(out)
    }) %>%
        bindEvent(input$load_UMI_data)

    working_output_filename <- reactive({
        if (input$UMAP_type == "multiple samples analysis (Seurat CCA)") {
            suffix <- paste0(
                input$multisample, "_",
                paste(input$plotting_sample, collapse = "&")
            )
        } else {
            suffix <- paste0(
                "Single_",
                paste(input$plotting_sample, collapse = "&")
            )
        }
        out <- paste0(
            format(Sys.time(), format = "%Y%m%d%H%M%S"), "_",
            input$plotting_feature, "_",
            suffix, "_",
            paste(input$color_range, collapse = "-"),
            ifelse(input$output_clear, "_Clear", ""),
            ".png"
        )
        return(out)
    })

    output$image <- renderImage({
        if (!dir.exists(output_folder)) {
            dir.create(output_folder, recursive = TRUE)
        }

        output_image_path <- file.path(output_folder, working_output_filename())
        plot_figure_out(
            plotting_df = get_plotting_df(),
            norm_factor_df = get_norm_factor_df(),
            UMI_Matrix = working_UMI_Matrix(),
            plotting_feature = input$plotting_feature,
            color_range = input$color_range,
            main = sub(
                pattern = ".png", replacement = "",
                substring(working_output_filename(), 16)
            ),
            output_clear = input$output_clear,
            output_path = output_image_path
        )

        # Return a list containing information about the image
        list(
            src = output_image_path,
            contentType = "image/png",
            width = 500
        )
    }) %>%
        bindEvent(input$figure_out)

    output$download_filename_setting <- renderUI(
        textInput(
            inputId = "download_filename",
            label = "Downloaded file name:",
            value = working_output_filename(),
            width = "500px"
        )
    )

    output$download_figure <-
        downloadHandler(
            filename = function() input$download_filename,
            content = function(file) {
                plot_figure_out(
                    plotting_df = get_plotting_df(),
                    norm_factor_df = get_norm_factor_df(),
                    UMI_Matrix = working_UMI_Matrix(),
                    plotting_feature = input$plotting_feature,
                    color_range = input$color_range,
                    main = sub(
                        pattern = ".png", replacement = "",
                        substring(working_output_filename(), 16)
                    ),
                    output_clear = input$output_clear,
                    output_path = file
                )
            }
        )
}
#