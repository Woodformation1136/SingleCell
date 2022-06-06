library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(Matrix)
library(scales)
library(slingshot)
library(forecast)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
# Output figure
output_png_figure <- function(plotting_function,
                              # x, y, col, main = "",
                              output_figure = FALSE,
                              output_path = "temp.png",
                              output_without_margin = FALSE,
                              ...) {
    if (output_figure) {
        png(output_path,
            pointsize = 12, res = 300,
            width = 32, height = 40, units = "cm"
        )
    }

    if (output_without_margin) {
        par(mai = c(0, 0, 0, 0))
    } else {
        par(mai = ori_par$mai)
    }

    plotting_function(
        output_figure = output_figure,
        output_path = output_path,
        output_without_margin = output_without_margin,
        ...
    )

    par(mai = ori_par$mai)

    if (output_figure) dev.off()
}

# Get pseudotime matrix along each lineage
get_pseudotime_matrix <- function(seurat_umap, sample_name) {
    plot_seurat_umap <- filter(seurat_umap, Sample == sample_name)
    n_cluster <- cluster_number_list[[sample_name]]

    if (sample_name == "TenX_Lch") {
        foo <- plot_seurat_umap %>%
            filter(CellrangerColor == "#9467BD") %>%
            arrange(UMAP_1)
        excluded_barcodes <- foo$SampleBarcode %>%
            head(floor(nrow(foo) * 5 / 100))
        plot_seurat_umap <- plot_seurat_umap %>%
            filter(!(SampleBarcode %in% excluded_barcodes))
    }

    if (sample_name == "TenX_PalChen2021") {
        foo <- plot_seurat_umap %>%
            filter(CellrangerColor == "#E377C2") %>%
            arrange(UMAP_1)
        excluded_barcodes <- foo$SampleBarcode %>%
            tail(floor(nrow(foo) * 5 / 100))
        plot_seurat_umap <- plot_seurat_umap %>%
            filter(!(SampleBarcode %in% excluded_barcodes))
    }

    if (sample_name == "MARSseq_Egr") {
        foo <- plot_seurat_umap %>%
            filter(CellrangerColor %in% c("#FF7F0F", "#F8E71C")) %>%
            arrange(UMAP_1)
        excluded_barcodes <- foo$SampleBarcode %>%
            tail(floor(nrow(foo) * 5 / 100))
        plot_seurat_umap <- plot_seurat_umap %>%
            filter(!(SampleBarcode %in% excluded_barcodes))
    }

    # Identify the psuedotime of each cell on each lineage
    lineage <- new("SlingshotDataSet")
    # str(lineage)

    ## Prepare the reducedDim
    rd <- as.matrix(plot_seurat_umap[, c("UMAP_1", "UMAP_2")])

    ## Prepare the clusterLabels
    cluster_factor <- factor(plot_seurat_umap$CellrangerColor)
    # plot(1:10,1:10,col=levels(cluster_factor),pch=20,cex=5)
    cl <- sapply(
        as.numeric(cluster_factor),
        function(i) {
            out <- rep(0, n_cluster)
            out[i] <- 1
            return(out)
        }
    ) %>% t()
    rownames(cl) <- rownames(rd)
    colnames(cl) <- seq(n_cluster)

    ## Prepare the lineages
    lin <- sapply(
        list(
            c("#FF7F0F", "#F8E71C", "#E377C2"),
            c("#9467BD", "#2AA02A", "#8C564B", "#1F77B4"),
            c("#9467BD", "#2AA02A", "#8C564B", "#D62728")
        ),
        function(lineage) {
            out <- match(lineage, levels(cluster_factor)) %>%
                as.character() %>%
                na.omit() %>%
                set_attributes(NULL)
            return(out)
        }
    ) %>%
        set_names(c("Lineage1", "Lineage2", "Lineage3"))

    ## Prepare the adjacency
    adc <- matrix(0, n_cluster, n_cluster)
    rownames(adc) <- seq(n_cluster)
    colnames(adc) <- seq(n_cluster)
    for (L in lin) {
        for (Ci in seq(length(L) - 1)) {
            adc[L[Ci], L[Ci + 1]] <- 1
            adc[L[Ci + 1], L[Ci]] <- 1
        }
    }

    ## Fill in the contents
    lineage@reducedDim <- rd
    lineage@clusterLabels <- cl
    lineage@lineages <- lin
    lineage@adjacency <- adc

    ## Plot the lineages
    # plot(plot_seurat_umap$UMAP_1,
    #      plot_seurat_umap$UMAP_2,
    #      pch = 20, cex = 0.3,
    #      col = plot_seurat_umap$CellrangerColor,
    #      xlab = "UMAP_1", ylab = "UMAP_2", main = "")
    # lines(lineage, lwd = 3, col = "black")

    ## Get lineage curves
    lineage <- getCurves(lineage, extend = "n", stretch = 0)

    pseudotime_matrix <- slingPseudotime(lineage)
    rownames(pseudotime_matrix) <-
        plot_seurat_umap %>%
        extract(, "SampleBarcode") %>%
        sub(pattern = paste0(sample_name, "_"), replacement = "")

    return(pseudotime_matrix)
}

# Get cellranger cluster color along each lineage
get_pseudotime_color <- function(seurat_umap, sample_name) {
    pseudotime_color <-
        seurat_umap %>%
        filter(Sample == sample_name) %>%
        extract(, "CellrangerColor")
    names(pseudotime_color) <-
        seurat_umap %>%
        filter(Sample == sample_name) %>%
        extract(, "SampleBarcode") %>%
        sub(pattern = paste0(sample_name, "_"), replacement = "")
    return(pseudotime_color)
}

# plot ortholog expression v.s. pseudotime
plot_log2_norm_ortho_UMI_on_pseudotime <- function(
    seurat_umap = seurat_umap,
    sample_name = "TenX_Ptr",
    log2_norm_ortho_UMI_matrix,
    ortholog_cluster_name  = paste0("Cluster_", plot_ortholog$OrthoGroupNo[4]),
    lineage_id = 1,
    ...
) {
    pseudotime_matrix <- get_pseudotime_matrix(seurat_umap, sample_name)
    pseudotime_color <- get_pseudotime_color(seurat_umap, sample_name)

    X <- pseudotime_matrix[, lineage_id] %>%
        na.omit() %>%
        sort()
    Y <- log2_norm_ortho_UMI_matrix[ortholog_cluster_name, names(X)]
    COL <- pseudotime_color[names(X)]

    plot(
        X, Y, col = COL, pch = 20, las = 1,
        main = "", xlab = "", ylab = "", ...
    )
    # main = paste(sample_name, "lineage", lineage_id),
    # xlab = "Pseudotime",
    # ylab = "Log2 normalized ortholog expression",
    # points(
    #     X, rep(par()$usr[3], length(X)),
    #     col = COL, pch = "|", cex = 0.5,
    #     xpd = TRUE
    # )
    # Y_ma <- ma(Y, ma_order)
    # Y_lowess <- lowess(X[!is.na(X)], Y[!is.na(X)], f = 1)
    # lines(X[!is.na(X)], Y[!is.na(X)])
}


# Set parameters ===============================================================
# Get input parameters from snakemake
input_integrated_rds <- snakemake@input$input_integrated_rds
input_dirs <- snakemake@input$input_dirs
log2_norm_ortho_dirs <- snakemake@input$log2_norm_ortho_dirs
plot_ortholog_csv <- snakemake@input$plot_ortholog_csv
cluster_number_list <- snakemake@params$cluster_number_list
color_list <- snakemake@params$color_list
output_dir <- snakemake@output$output_dir

# plot_ortholog_csv <- "rawdata/20220307_orthogroup_list.csv"
# input_integrated_rds <- "results/Integrated_rds/Species4_Sample5.rds"
# log2_norm_ortho_dirs <- "results/Ortholog_UMI_matrix/Each_log2_norm"
# output_dir <- "results/Integrated_figures/UMAP_by_cellranger_cluster_on_sample/Species4_Sample5"
# cluster_number_list <- list(
#     TenX_Ptr = 10,
#     TenX_Lch = 10,
#     TenX_PalChen2021 = 20,
#     MARSseq_Egr = 10,
#     MARSseq_Tar = 18
# )
# color_list <- list(
#     TenX_Ptr = c(
#         "#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C",
#         "#9467BD", "#D62728", "#E377C2", "#9B9B9B", "#4B4B4B"
#     ),
#     TenX_Lch = c(
#         "#1F77B4", "#8C564B", "#D62728", "#FF7F0F", "#F8E71C",
#         "#E377C2", "#9467BD", "#2AA02A", "#63EE9B", "#9B9B9B"
#     ),
#     TenX_PalChen2021 = c(
#         "#1F77B4", "#8C564B", "#E377C2", "#1F77B4", "#F8E71C",
#         "#D62728", "#D62728", "#2AA02A", "#1F77B4", "#FF7F0F",
#         "#E377C2", "#2AA02A", "#9B9B9B", "#9467BD", "#2AA02A",
#         "#D62728", "#E377C2", "#FF7F0F", "#4B4B4B", "#8C564B"
#     ),
#     MARSseq_Egr = c(
#         "#9467BD", "#F8E71C", "#8C564B", "#1F77B4", "#2AA02A",
#         "#FF7F0F", "#D62728", "#E377C2", "#4B4B4B", "#9B9B9B"
#     ),
#     MARSseq_Tar = c(
#         "#1F77B4", "#E377C2", "#55A3FF", "#8C564B", "#9467BD",
#         "#FF7F0F", "#46CBE5", "#F8E71C", "#2AA02A", "#9B9B9B",
#         "#9B9B9B", "#9B9B9B", "#9B9B9B", "#9B9B9B", "#9B9B9B",
#         "#9B9B9B", "#9B9B9B", "#9B9B9B"
#     )
# )
# input_dirs <- c(
#     "results/SCseq_TenX_Ptr", "results/SCseq_TenX_PalChen2021",
#     "results/SCseq_MARSseq_Egr", "results/SCseq_MARSseq_Tar",
#     "results/SCseq_TenX_Lch"
# )

input_samples <- sapply(
    input_dirs,
    function(x) {
        suffix <- sub(pattern = "results/SCseq_", replacement = "", x)
        return(suffix)
    }
)
n_samples <- length(input_samples)
input_clustering_table_path <- sapply(
    input_samples,
    function(x) {
        out <- paste0(
            "results/SCseq_", x, "/cellranger_reanalysis_", x,
            "/outs/analysis/clustering/kmeans_", cluster_number_list[[x]],
            "_clusters/clusters.csv"
        )
        return(out)
    }
)
input_clustering_table_info <- cbind(input_samples, input_clustering_table_path)


# Implement program ============================================================
# Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input cellranger clustering data
read_clustering_table <- function(table_info) {
    sample_name <- table_info[1]
    path <- table_info[2]
    df <- read.csv(path)
    df$SampleBarcode <- paste(sample_name, df$Barcode, sep = "_")
    df$CellrangerCluster <- df$Cluster
    return(df[, c("SampleBarcode", "CellrangerCluster")])
}
cellranger_cluster_table <- do.call(
    rbind,
    apply(input_clustering_table_info, 1, read_clustering_table)
)
rownames(cellranger_cluster_table) <- NULL

# Input integrated data
combined_object <- readRDS(input_integrated_rds)
seurat_umap <- as.data.frame(
    combined_object@reductions$umap@cell.embeddings
)
seurat_umap$SampleBarcode <- rownames(seurat_umap)
seurat_umap$Sample <-
    seurat_umap$SampleBarcode %>%
    strsplit(split = "_") %>%
    sapply(extract, simplify = FALSE, 1:2) %>%
    sapply(paste, collapse = "_")

# Merge seurat umap data and cellranger clustering data
seurat_umap <- left_join(
    seurat_umap,
    cellranger_cluster_table,
    by = "SampleBarcode"
)
stopifnot(
    all(!is.na(seurat_umap$CellrangerCluster))
)

# Get cellranger color from color list
seurat_umap$CellrangerColor <- sapply(
    seq(nrow(seurat_umap)),
    function(i) {
        sample_name <- seurat_umap$Sample[i]
        cellranger_cluster <- seurat_umap$CellrangerCluster[i]
        out <- color_list[[sample_name]][cellranger_cluster]
        return(out)
    }
)

# Input ortholog table for plotting
plot_ortholog <- read.csv(plot_ortholog_csv)

# Wrap function for plotting
plotting_function <- function(
    ortholog_cluster_id,
    ortholog_common_name,
    ...
) {
    par(mar = c(0, 0, 0, 0) + 0.1)

    n_filler <- n_samples * 5
    layout_matrix <- rbind(
        matrix(
            c(
                1, 1, 1, 1, 1,
                2, 2, 3, 4, 5
            ),
            ncol = 5, byrow = TRUE
        ),
        cbind(
            matrix(
                5 + seq(n_samples * 2),
                ncol = 2, byrow = FALSE
            ),
            matrix(
                5 + n_samples * 2 + seq(n_samples * 3),
                ncol = 3, byrow = TRUE
            )
        ),
        c(1, 1, 2, 3, 4) + 5 + n_samples * 5
    )
    layout(
        layout_matrix,
        heights = c(0.2, 0.2, rep(1, n_samples), 0.2),
        widths = c(0.2, 0.2, 1, 1, 1)
    )

    plot.new()
    text(0, 0.5, ortholog_common_name, cex = 4, adj = c(0, 0.5))

    plot.new()

    lineage_name <- c(
        "Ray lineage",
        "Fusiform lineage_vessel",
        "Fusiform lineage_fiber"
    )
    for (i in seq_along(lineage_name)) {
        plot.new()
        text(0.5, 0.5, lineage_name[i], cex = 2.5, adj = c(0.5, 0.5))
    }

    for (i in seq_along(input_samples)) {
        plot.new()
        # For convenience of modification...
        # text(
        #     0.5, 0.5, input_samples[i],
        #     cex = 2.5, adj = c(0.5, 0.5), srt = 90
        # )
    }

    for (i in seq_along(input_samples)) {
        plot.new()
        text(
            0.5, 0.5, "log2 normalized UMI",
            cex = 2, adj = c(0.5, 0.5), srt = 90
        )
    }

    par(mar = c(2, 2, 2.5, 1) + 0.1)
    for (j in seq(n_samples)) {
        for (i in seq(3)) {
            log2_norm_ortho_UMI_matrix <-
                readRDS(
                    paste0(log2_norm_ortho_dirs, "/", input_samples[j], ".rds")
                )
            ortholog_cluster_name <- paste0("Cluster_", ortholog_cluster_id)
            all_ortholog_cluster_name <- rownames(log2_norm_ortho_UMI_matrix)
            if (
                ortholog_cluster_name %in% all_ortholog_cluster_name &
                !(i == 3 & regexpr("Tar", input_samples[j]) != -1)
            ) {
                ylim_max <-
                    max(log2_norm_ortho_UMI_matrix[ortholog_cluster_name, ])
                plot_log2_norm_ortho_UMI_on_pseudotime(
                    seurat_umap = seurat_umap,
                    sample_name = input_samples[j],
                    log2_norm_ortho_UMI_matrix = log2_norm_ortho_UMI_matrix,
                    ortholog_cluster_name = ortholog_cluster_name,
                    lineage_id = i,
                    ylim = c(0, ylim_max),
                    axes = FALSE
                )

                par_xaxp <- par()$xaxp
                par_xaxp_gap <- (par_xaxp[2] - par_xaxp[1]) / par_xaxp[3]
                par_xaxp_new <- c(
                    par_xaxp[1] - par_xaxp_gap,
                    par_xaxp[2] + par_xaxp_gap,
                    par_xaxp[3] + 2
                )
                axis(1, cex.axis = 1.5, lwd = 2, xaxp = par_xaxp_new)

                if (i == 1) {
                    axis(2, cex.axis = 1.5, lwd = 2, las = 1)
                    box(bty = "L", lwd = 2)
                }
            } else {
                plot.new()
            }
        }
    }

    par(mar = c(0, 0, 0, 0) + 0.1)
    plot.new()

    for (i in 1:3) {
        plot.new()
        text(0.5, 0.5, "Pseudotime", cex = 2, adj = c(0.5, 0.5))
    }

    par(ori_par)
}

# Output figure
for (i_ortholog in seq(nrow(plot_ortholog))) {
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path = paste0(
            output_dir, "/Ortholog_",
            plot_ortholog$OrthoGroupNo[i_ortholog],
            "_",
            gsub("/", " ", plot_ortholog$CommonName[i_ortholog]),
            ".png"
        ),
        output_without_margin = FALSE,
        ortholog_cluster_id = plot_ortholog$OrthoGroupNo[i_ortholog],
        ortholog_common_name = plot_ortholog$CommonName[i_ortholog]
    )
}
