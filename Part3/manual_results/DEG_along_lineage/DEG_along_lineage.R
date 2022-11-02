# devtools::install_github("jokergoo/ComplexHeatmap")
library(slingshot)
library(magrittr)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
oriPar <- par(no.readonly = TRUE)


# Define functions =============================================================
# Plot heatmap of expression along each lineage
heatmap_along_lineage <- function(
    plot_df = plotting_TenX_Ptr_df,
    plot_features_list = lineage_gene_list,
    plot_UMI_matrix = t(plot_norm_UMI_matrix),
    plot_pseudotime_list = Ptr_lineage_list,
    plot_filename,
    k_nn_color = 10,
    bandwidth = 0.05,
    steplength = 0.005,
    output_figure = TRUE,
    upper_col = "#CE5324",
    lower_col = "#f3d3b8",
    bar_col = "#FDC441"
) {
    plot_features <- unlist(plot_features_list)
    plot_UMI_matrix <- plot_UMI_matrix[plot_features, ]

    plot_info_list <- sapply(
        simplify = FALSE,
        plot_pseudotime_list,
        function(pseudotime) {
            barcode_vector <- names(pseudotime)
            color_vector <- plot_df[barcode_vector, ]$Color
            UMI_matrix <- plot_UMI_matrix[, barcode_vector]
            sampling_pseudotime <-
                seq(min(pseudotime), max(pseudotime), steplength)
            estimated_exp_matrix <- apply(
                UMI_matrix, 1,
                function(exp_vector) {
                    out <- ksmooth(
                        pseudotime,
                        exp_vector,
                        kernel = "normal",
                        bandwidth = bandwidth,
                        x.points = sampling_pseudotime
                    )$y
                    return(log2(out + 1))
                }
            )
            estimated_color_vector <- sapply(
                sampling_pseudotime,
                function(x) {
                    tmp <- 1
                    repeat {
                        is_nn <-
                            (pseudotime >= x - bandwidth * tmp) &
                            (pseudotime <= x + bandwidth * tmp)
                        if (sum(is_nn) >= k_nn_color) break
                        tmp <- tmp + 1
                    }
                    nncolor <- color_vector[is_nn]
                    out <- names(which.max(table(nncolor)))
                    return(out)
                }
            )
            return(list(
                sampling_pseudotime = sampling_pseudotime,
                estimated_exp_matrix = estimated_exp_matrix,
                estimated_color_vector = estimated_color_vector
            ))
        }
    )


    plot_expression_matrix <-
        do.call(
            rbind,
            sapply(plot_info_list, extract2, "estimated_exp_matrix")
        ) %>%
        t()
    plot_expression_matrix <- rbind(plot_expression_matrix, 0)
    plot_expression_matrix <- rbind(plot_expression_matrix, 0)
    colnames(plot_expression_matrix) <-
        paste0("Cell_", seq_len(ncol(plot_expression_matrix)))
    rownames(plot_expression_matrix) <-
        paste0("Feature_", seq_len(nrow(plot_expression_matrix)))

    feature_group <- rep(
        names(plot_features_list),
        times = sapply(plot_features_list, length)
    )
    feature_group <- c(feature_group, "(= _ =)", "(= _ =)")
    plot_mean_UMI <- apply(plot_UMI_matrix, 1, mean)
    anno_row <- data.frame(
        Group = feature_group,
        row.names = rownames(plot_expression_matrix)
    )

    EachMeanUMI_df <- sapply(
        simplify = FALSE,
        plot_pseudotime_list,
        function(x) rowMeans(plot_UMI_matrix[, names(x)]) / plot_mean_UMI
    ) %>%
    rev() %>%
    as.data.frame()
    EachMeanUMI_df <-
        rbind(
            EachMeanUMI_df,
            rep(max(EachMeanUMI_df), ncol(EachMeanUMI_df)),
            rep(min(EachMeanUMI_df), ncol(EachMeanUMI_df))
        )
    anno_row <- cbind(anno_row, EachMeanUMI_df)

    anno_col <- data.frame(
        Cluster = unlist(
            sapply(plot_info_list, extract2, "estimated_color_vector")
        ),
        row.names = colnames(plot_expression_matrix)
    )
    ann_colors <- list(
        Group = set_names(
            c("#26C3F7", "#F969A8", "#1B9E77", "#DDDDDD"),
            c(names(plot_features_list), "(= _ =)")
        ),
        Cluster = set_names(
            unique(anno_col$Cluster),
            unique(anno_col$Cluster)
        )
    )
    ann_colors <-
        c(
            ann_colors,
            sapply(
                simplify = FALSE, USE.NAMES = TRUE,
                colnames(EachMeanUMI_df),
                function(tmp) c(lower_col, upper_col)
            )
        )

    scale_plot_expression_matrix <-
        plot_expression_matrix %>%
        apply(1, function(x) {
            out <- punif(x, min(x, na.rm = TRUE), max(x, na.rm = TRUE))
            out <- pmax(pmin(out * 1.4 - 0.2, 1), 0) * 100
            return(out)
        }) %>%
        t() %>%
        set_colnames(colnames(plot_expression_matrix))

    n_sampling_pseudotime <- sapply(
        plot_info_list, function(x) length(x[["sampling_pseudotime"]])
    )
    out_heatmap <-
        pheatmap(
            scale_plot_expression_matrix,
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            color = cividis(100),
            breaks = seq(0, max(scale_plot_expression_matrix, na.rm = TRUE),
                        length.out = 100),
            border_color = NA,
            gaps_row = cumsum(sapply(plot_features_list, length)),
            gaps_col = cumsum(n_sampling_pseudotime),
            show_rownames = FALSE,
            show_colnames = FALSE,
            fontsize_row = 8,
            annotation_row = anno_row,
            annotation_col = anno_col,
            annotation_colors = ann_colors,
            na_col = "gray80"
        )

    out_heatmap <-
        out_heatmap +
        HeatmapAnnotation(
            which = "row",
            gap = unit(10, "points"),
            MF = anno_barplot(
                EachMeanUMI_df$Vessel,
                gp = gpar(col = bar_col, fill = bar_col)
            ),
            MR = anno_barplot(
                EachMeanUMI_df$Ray,
                gp = gpar(col = bar_col, fill = bar_col)
            )
        )

    if (output_figure) {
        png(plot_filename, width = 40, height = 20, units = "cm", res = 600)
        draw(out_heatmap)
        dev.off()
    }
}


# Set parameters ===============================================================
plotting_TenX_Ptr_csv <- file.path()
norm_factor_TenX_Ptr_csv <- file.path()
geneUMI_TenX_Ptr_csv <- file.path()

lineage_gene_list <- list(
    "Fusiform lineage" = c(
        "Potri.001G240900.v4.1",
        "Potri.001G102400.v4.1",
        "Potri.004G172000.v4.1",
        "Potri.001G186700.v4.1",
        "Potri.009G092300.v4.1",
        "Potri.004G167700.v4.1",
        "Potri.001G040200.v4.1",
        "Potri.011G057500.v4.1",
        "Potri.012G110750.v4.1",
        "Potri.009G107050.v4.1",
        "Potri.012G047800.v4.1",
        "Potri.016G088700.v4.1",
        "Potri.011G134900.v4.1",
        "Potri.002G065200.v4.1",
        "Potri.005G101800.v4.1",
        "Potri.005G228600.v4.1"
    ),
    "Ray cell lineage" = c(
        "Potri.T011200.v4.1",
        "Potri.008G144500.v4.1",
        "Potri.004G196000.v4.1",
        "Potri.003G015200.v4.1",
        "Potri.002G223100.v4.1",
        "Potri.001G218800.v4.1",
        "Potri.004G108320.v4.1",
        "Potri.005G195600.v4.1"
    ),
    "MYB" = c(
        "Potri.009G134000.v4.1",
        "Potri.004G174400.v4.1"
    )
)


# Implement analyses ===========================================================
plotting_TenX_Ptr_df <- read.csv(plotting_TenX_Ptr_csv, row.names = 1)
norm_factor_TenX_Ptr_df <- read.csv(norm_factor_TenX_Ptr_csv, row.names = 1)
geneUMI_TenX_Ptr_df <- read.csv(geneUMI_TenX_Ptr_csv, row.names = 1)

norm_factor_TenX_Ptr_df <-
    norm_factor_TenX_Ptr_df[rownames(plotting_TenX_Ptr_df), "norm_factor"]
geneUMI_TenX_Ptr_df <-
    geneUMI_TenX_Ptr_df[rownames(plotting_TenX_Ptr_df), ]
plot_norm_UMI_matrix <-
    apply(geneUMI_TenX_Ptr_df, 2, multiply_by, norm_factor_TenX_Ptr_df)

#Identify the pseudotime of each cell on each lineage
Ptr_lineage <- new("SlingshotDataSet")

##Prepare the reducedDim
rd <- as.matrix(plotting_TenX_Ptr_df[, c("UMAP.1", "UMAP.2")])

##Prepare the clusterLabels
cluster_factor <- factor(
    plotting_TenX_Ptr_df$Color,
    levels = c(
        "#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C",
        "#9467BD", "#D62728", "#E377C2", "#9B9B9B", "#4B4B4B"
    )
)

cl <- sapply(
    as.numeric(cluster_factor),
    function(i) {
        out <- rep(0, 10)
        out[i] <- 1
        return(out)
    }
) %>% t()
rownames(cl) <- rownames(plotting_TenX_Ptr_df)
colnames(cl) <- 1:10

##Prepare the lineages
lin <- list(
    Lineage1 = c("3", "5", "8"),
    Lineage2 = c("6", "4", "2", "1"),
    Lineage3 = c("6", "4", "2", "7")
)

##Prepare the adjacency
adc <- matrix(0, 10, 10)
rownames(adc) <- 1:10
colnames(adc) <- 1:10
adc[3, 5] <- adc[5, 3] <- 1
adc[5, 8] <- adc[8, 5] <- 1
adc[6, 4] <- adc[4, 6] <- 1
adc[4, 2] <- adc[2, 4] <- 1
adc[2, 1] <- adc[1, 2] <- 1
adc[2, 7] <- adc[7, 2] <- 1

##Fill in the contents
Ptr_lineage@reducedDim <- rd
Ptr_lineage@clusterLabels <- cl
Ptr_lineage@lineages <- lin
Ptr_lineage@adjacency <- adc

##Get lineage curves
Ptr_lineage <- getCurves(Ptr_lineage, extend = "y")

##Get pseudotime matrix
Ptr_pseudotime_matrix <- slingPseudotime(Ptr_lineage)
Ptr_lineage_list <-
    sapply(
        1:3,
        function(lid) {
            lineage_barcodes <-
                Ptr_pseudotime_matrix[, lid] %>%
                extract(!is.na(.)) %>%
                sort()
            return(lineage_barcodes)
        }
    )
names(Ptr_lineage_list) <- c("Ray", "Vessel", "Fiber")


heatmap_along_lineage(
    plot_df = plotting_TenX_Ptr_df,
    plot_features_list = lineage_gene_list,
    plot_UMI_matrix = t(plot_norm_UMI_matrix),
    plot_pseudotime_list = Ptr_lineage_list,
    plot_filename = file.path(),
    k_nn_color = 10,
    bandwidth = 0.35,
    steplength = 0.05,
    output_figure = TRUE
)

heatmap_along_lineage(
    plot_df = plotting_TenX_Ptr_df,
    plot_features_list = lineage_gene_list,
    plot_UMI_matrix = t(plot_norm_UMI_matrix),
    plot_pseudotime_list = Ptr_lineage_list,
    plot_filename = file.path(),
    k_nn_color = 10,
    bandwidth = 0.35,
    steplength = 0.05,
    output_figure = TRUE,
    upper_col = "#D25030",
    lower_col = "#E8AA98",
    bar_col = "#CB6543"
)
heatmap_along_lineage(
    plot_df = plotting_TenX_Ptr_df,
    plot_features_list = lineage_gene_list,
    plot_UMI_matrix = t(plot_norm_UMI_matrix),
    plot_pseudotime_list = Ptr_lineage_list,
    plot_filename = file.path(),
    k_nn_color = 10,
    bandwidth = 0.35,
    steplength = 0.05,
    output_figure = TRUE,
    upper_col = "#34967F",
    lower_col = "#AEDCD2",
    bar_col = "#42a991"
)
heatmap_along_lineage(
    plot_df = plotting_TenX_Ptr_df,
    plot_features_list = lineage_gene_list,
    plot_UMI_matrix = t(plot_norm_UMI_matrix),
    plot_pseudotime_list = Ptr_lineage_list,
    plot_filename = file.path(),
    k_nn_color = 10,
    bandwidth = 0.35,
    steplength = 0.05,
    output_figure = TRUE,
    upper_col = "#348028",
    lower_col = "#bfe689",
    bar_col = "#68a033"
)
heatmap_along_lineage(
    plot_df = plotting_TenX_Ptr_df,
    plot_features_list = lineage_gene_list,
    plot_UMI_matrix = t(plot_norm_UMI_matrix),
    plot_pseudotime_list = Ptr_lineage_list,
    plot_filename = file.path(),
    k_nn_color = 10,
    bandwidth = 0.35,
    steplength = 0.05,
    output_figure = TRUE,
    upper_col = "#276696",
    lower_col = "#8BCED8",
    bar_col = "#3392B6"
)
