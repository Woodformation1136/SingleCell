library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(Matrix)
library(scales)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
# Output figure
output_png_figure <- function(
    plotting_function,
    # x, y, col, main = "",
    output_figure = FALSE,
    output_path = "temp.png",
    output_without_margin = FALSE,
    ...
) {
    if (output_figure) {
        png(output_path,
            pointsize = 10, res = 300,
            width = 20, height = 15, units = "cm")
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

# Plot seurat UMAP colored on each seurat cluster
plot_on_seurat_cluster <- function(
    seurat_umap,
    seurat_cluster_info,
    use_cellranger_cluster_color,
    output_without_margin,
    output_without_number,
    ...
) {
    x <- seurat_umap$UMAP_1
    y <- seurat_umap$UMAP_2
    if (use_cellranger_cluster_color) {
        col <- with(
            seurat_cluster_info,
            color[match(seurat_umap$SeuratCluster, SeuratCluster)]
        )
    } else {
        n_col <- length(unique(seurat_umap$SeuratCluster))
        col_vector <- scales::hue_pal()(n_col)
        col <- col_vector[seurat_umap$SeuratCluster]
    }

    randam_order <- sample(length(x))
    plot(
        x[randam_order], y[randam_order],
        col = col[randam_order], pch = 20, cex = 0.3,
        xlab = "UMAP_1", ylab = "UMAP_2",
        main = ifelse(output_without_margin, "", "main"),
        axes = !output_without_margin, las = 1
    )
    # for (i in seq(nrow(seurat_cluster_info))) {
    #     text(
    #         seurat_cluster_info$UMAP_1[i],
    #         seurat_cluster_info$UMAP_2[i],
    #         seurat_cluster_info$SeuratCluster[i],
    #         cex = 1.5
    #     )
    # }
    if (!output_without_number) {
        text(
            seurat_cluster_info$UMAP_1,
            seurat_cluster_info$UMAP_2,
            seurat_cluster_info$SeuratCluster,
            cex = 2.5
        )
    }
}

# Wrap function of plotting seurat UMAP colored on each sample
plotting_function <- function(...) {
    plot_on_seurat_cluster(
        seurat_umap = seurat_umap,
        seurat_cluster_info = seurat_cluster_info,
        ...
    )
}

# Set parameters ===============================================================
# Get input parameters from command line
input_integrated_rds <- snakemake@input$input_integrated_rds
input_dirs <- snakemake@input$input_dirs
cluster_number_list <- snakemake@params$cluster_number_list
color_list <- snakemake@params$color_list
output_dir <- snakemake@output$output_dir

input_samples <- sapply(
    X = input_dirs,
    FUN = sub,
    pattern = "results/SCseq_",
    replacement = ""
)
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


# Implementation ===============================================================
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
    df$CellrangerColor <- color_list[[sample_name]][df$CellrangerCluster]

    return(df[, c("SampleBarcode", "CellrangerCluster", "CellrangerColor")])
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
seurat_umap$SeuratCluster <- as.numeric(
    combined_object@meta.data$seurat_clusters
)
seurat_umap$SampleBarcode <- rownames(seurat_umap)
seurat_umap$Sample <-
    seurat_umap$SampleBarcode %>%
    strsplit(split = "_") %>%
    sapply(extract, simplify = FALSE, 1:2) %>%
    sapply(paste, collapse = "_")
n_samples <- length(unique(seurat_umap$Sample))

# Merge seurat umap data and cellranger clustering data
seurat_umap <- left_join(
    seurat_umap,
    cellranger_cluster_table,
    by = "SampleBarcode"
)
stopifnot(
    all(!is.na(seurat_umap$CellrangerCluster))
)

# Construct color and center coordinates for each seurat cluster
seurat_cluster_info <- seurat_umap %>%
    group_by(SeuratCluster) %>%
    summarise(
        color = names(which.max(table(CellrangerColor))),
        UMAP_1 = mean(UMAP_1),
        UMAP_2 = mean(UMAP_2)
    )

## Output figure
output_png_figure(
    plotting_function,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/BySeuratCluster_Elected_color.png"
        ),
    output_without_margin = FALSE,
    use_cellranger_cluster_color = TRUE,
    output_without_number = TRUE
)
output_png_figure(
    plotting_function,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/BySeuratCluster_Seperated_color.png"
        ),
    output_without_margin = FALSE,
    use_cellranger_cluster_color = FALSE,
    output_without_number = FALSE
)
output_png_figure(
    plotting_function,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/BySeuratCluster_Elected_color_Clear.png"
        ),
    output_without_margin = TRUE,
    use_cellranger_cluster_color = TRUE,
    output_without_number = TRUE
)
output_png_figure(
    plotting_function,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/BySeuratCluster_Seperated_color_Clear.png"
        ),
    output_without_margin = TRUE,
    use_cellranger_cluster_color = FALSE,
    output_without_number = FALSE
)