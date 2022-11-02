library(Seurat)
library(magrittr)
library(reticulate)


# Define functions =============================================================


# Set parameters ===============================================================
# Get input parameters from command line
input_integration_rds <- snakemake@input$integration_rds
output_integrated_umap_csv <- snakemake@output$integrated_umap_csv
seed.use <- as.numeric(snakemake@wildcards$seed_use)
min.dist <- as.numeric(snakemake@wildcards$min_dist)
n.neighbors <- as.numeric(snakemake@wildcards$n_neighbors)


# Implementation ===============================================================
# Input combined object
combined_object <- readRDS(input_integration_rds)


# Run UMAP using python module umap (umap-learning in Seurat)
umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
sklearn <- reticulate::import("sklearn", delay_load = TRUE)
random.state <-
    sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
umap.args <- list(
    n_neighbors = as.integer(n.neighbors),
    n_components = 2L, metric = "correlation",
    n_epochs = NULL, learning_rate = 1,
    min_dist = min.dist, spread = 1, set_op_mix_ratio = 1,
    local_connectivity = 1L, repulsion_strength = 1,
    negative_sample_rate = 5, random_state = random.state,
    a = NULL, b = NULL, metric_kwds = NULL, angular_rp_forest = FALSE,
    verbose = TRUE, init = "random"
    # transform_seed = seed.use, init = "random"
)
umap <- do.call(what = umap_import$UMAP, args = umap.args)
base_matrix <- Embeddings(combined_object[["pca"]])[, 1:30]
umap_df <- cbind(
    Barcode = rownames(base_matrix),
    umap$fit_transform(base_matrix) %>%
    data.frame() %>%
    set_colnames(c("UMAP.1", "UMAP.2"))
)
write.csv(
    umap_df,
    file = output_integrated_umap_csv,
    row.names = FALSE,
    quote = FALSE
)
