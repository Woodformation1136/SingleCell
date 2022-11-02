suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(Matrix)
})

# Define functions =============================================================


# Set parameters ===============================================================
# Get input parameters from command line
input_integrated_rds <- snakemake@input$integrated_rds
input_integrated_umap_csv <- snakemake@input$integrated_umap_csv
output_plotting_csv <- snakemake@output$plotting_csv
sample_name <- snakemake@params$sample
species <- snakemake@params$species


# Implementation ===============================================================
# Create output directory
if (!dir.exists(dirname(output_plotting_csv))) {
    dir.create(dirname(output_plotting_csv), recursive = TRUE)
}

# Input integrated data
combined_object <- readRDS(input_integrated_rds)

# Get information from combined_object and store as data.frame
plotting_df <- read.csv(input_integrated_umap_csv)
feature_sample_barcode <- strsplit(plotting_df$Barcode, split = "_")

if (any(sapply(feature_sample_barcode, length) != 4)) {
    stop(
        "Not all barcode names from integrated result follow ",
        "the pattern 'feature_platform_batch_barcode', ",
        "so get_seuratCCA_results_into_csv.R need to been revised."
    )
}
feature_type <- feature_sample_barcode[[1]][1]
plotting_df$Sample <-
    feature_sample_barcode %>%
    sapply(extract, simplify = FALSE, 2:3) %>%
    sapply(paste, collapse = "_")
plotting_df$Barcode <-
    feature_sample_barcode %>%
    sapply(extract, simplify = TRUE, 4)
plotting_df$Cluster <-
    as.numeric(as.character(combined_object@meta.data$seurat_clusters))
plotting_df$Species <-
    species[match(plotting_df$Sample, sample_name)]

# Reorganize plotting_df
new_col_order <-
    c("Barcode", "Sample", "Species", "UMAP.1", "UMAP.2", "Cluster")
plotting_df <- plotting_df[, new_col_order]
rownames(plotting_df) <- NULL

write.csv(
    plotting_df,
    file = output_plotting_csv,
    row.names = FALSE, quote = FALSE
)