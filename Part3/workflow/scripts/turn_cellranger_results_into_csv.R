suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(Matrix)
})

# Define functions =============================================================
# Read UMI matrix from 10X format
read_UMI_matrix <- function(matrix_dir) {
    barcode_path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features_path <- paste0(matrix_dir, "features.tsv.gz")
    matrix_path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- t(readMM(file = matrix_path))
    feature_names <- read.delim(
        features_path,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    barcode_names <- read.delim(
        barcode_path,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(mat) <- barcode_names[, 1]
    colnames(mat) <- feature_names[, 2]
    return(mat)
}


# Set parameters ===============================================================
# Get input parameters from snakemake
input_folder <- snakemake@input$folder
total_UMI_count_cutoff_on_barcode <-
    snakemake@params$total_UMI_count_cutoff_on_barcode
n_cluster <- snakemake@params$n_cluster
cluster_color <- snakemake@params$cluster_color
output_geneUMI_csv <- snakemake@output$geneUMI_csv
output_plotting_csv <- snakemake@output$plotting_csv
output_norm_factor_csv <- snakemake@output$norm_factor_csv


# Implementation ===============================================================
# Create output directory
if (!dir.exists(dirname(output_geneUMI_csv))) {
    dir.create(dirname(output_geneUMI_csv), recursive = TRUE)
}

# Get UMI count data -----------------------------------------------------------
UMI_matrix_dir <- paste0(input_folder, "/outs/filtered_feature_bc_matrix/")
UMI_matrix <- read_UMI_matrix(UMI_matrix_dir)

# Keep barcodes with total UMI count >= total_UMI_count_cutoff_on_barcode
total_UMI_count <- rowSums(UMI_matrix)
filtered_UMI_matrix <-
    UMI_matrix[total_UMI_count >= total_UMI_count_cutoff_on_barcode, ]

# Convert into data.frame
UMI_df <- cbind(
    Barcode = rownames(filtered_UMI_matrix),
    data.frame(as.matrix(filtered_UMI_matrix))
)

write.csv(
    UMI_df,
    file = output_geneUMI_csv,
    row.names = FALSE,
    quote = FALSE
)


# Get plotting data ------------------------------------------------------------
cluster_csv_path <- paste0(
    input_folder,
    "/outs/analysis/clustering/kmeans_",
    n_cluster,
    "_clusters/clusters.csv"
)
cluster_df <- read.csv(cluster_csv_path)

umap_csv_path <- paste0(
    input_folder,
    "/outs/analysis/umap/2_components/projection.csv"
)
umap_df <- read.csv(umap_csv_path)

plotting_df <- full_join(cluster_df, umap_df, by = "Barcode")
plotting_df$Color <- cluster_color[plotting_df$Cluster]

# Keep barcodes with total UMI count >= total_UMI_count_cutoff_on_barcode
plotting_df <- plotting_df[match(UMI_df$Barcode, plotting_df$Barcode), ]

write.csv(
    plotting_df, file = output_plotting_csv,
    row.names = FALSE, quote = FALSE
)


# Get normalization factors ----------------------------------------------------
norm_factors_df <- data.frame(
    Barcode = rownames(filtered_UMI_matrix),
    norm_factor = 10000 / rowSums(filtered_UMI_matrix)
)

write.csv(
    norm_factors_df,
    file = output_norm_factor_csv,
    row.names = FALSE, quote = FALSE
)
