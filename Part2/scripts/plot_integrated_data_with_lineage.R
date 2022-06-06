library(RColorBrewer)
library(magrittr)
library(dplyr)
library(Matrix)
library(scales)
library(slingshot)
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

# Get lineage along cellranger clusters on seurat UMAP
get_lineage_object <- function(seurat_umap, sample_name) {
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

    return(lineage)
}

# Plot seurat UMAP with lineage
plot_lineage_on_Seurat_UMAP <- function(seurat_umap,
                                        sample_name,
                                        lineage_object,
                                        lineage_id = 1,
                                        ...) {
    plot_seurat_umap <- filter(seurat_umap, Sample == sample_name)
    plot(plot_seurat_umap$UMAP_1,
        plot_seurat_umap$UMAP_2,
        pch = 20, cex = 0.3,
        col = plot_seurat_umap$CellrangerColor,
        xlab = "UMAP_1", ylab = "UMAP_2", main = ""
    )
    lines(
        lineage_object@metadata$curves[[lineage_id]],
        lwd = 3, col = "black"
    )
}


# Set parameters ===============================================================
# Get input parameters from snakemake
input_integrated_rds <- snakemake@input$input_integrated_rds
input_dirs <- snakemake@input$input_dirs
cluster_number_list <- snakemake@params$cluster_number_list
color_list <- snakemake@params$color_list
output_dir <- snakemake@output$output_dir

# input_integrated_rds <- "results/Integrated_rds/Species4_Sample5.rds"
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
#     "results/SCseq_TenX_Lch", "results/SCseq_MARSseq_Egr",
#     "results/SCseq_MARSseq_Tar"
# )

input_samples <- sapply(
    input_dirs,
    function(x) {
        suffix <- sub(pattern = "results/SCseq_", replacement = "", x)
        return(suffix)
    }
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
sample_vector <- unique(seurat_umap$Sample)
n_samples <- length(sample_vector)

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

# Wrap function of plotting seurat UMAP with lineage
plotting_function <- function(
    seurat_umap,
    sample_name,
    ...
) {
    lineage_object <- get_lineage_object(
        seurat_umap = seurat_umap,
        sample_name = sample_name
    )
    plot_lineage_on_Seurat_UMAP(
        seurat_umap = seurat_umap,
        sample_name = sample_name,
        lineage_object = lineage_object,
        ...
    )
}

# Output figure
for (sample_name in sample_vector) {
    for (lineage_id in 1:3) {
        output_png_figure(
            plotting_function,
            output_figure = TRUE,
            output_path = paste0(
                output_dir, "/ByCellrangerCluster_", sample_name,
                "_lineage", lineage_id, ".png"
            ),
            output_without_margin = FALSE,
            seurat_umap = seurat_umap,
            sample_name = sample_name,
            lineage_id = lineage_id
        )
        output_png_figure(
            plotting_function,
            output_figure = TRUE,
            output_path = paste0(
                output_dir, "/ByCellrangerCluster_", sample_name,
                "_lineage", lineage_id, "_Clear.png"
            ),
            output_without_margin = TRUE,
            seurat_umap = seurat_umap,
            sample_name = sample_name,
            lineage_id = lineage_id
        )
    }
}
