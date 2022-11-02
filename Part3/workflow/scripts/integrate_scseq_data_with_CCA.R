library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(tools)


# Define functions =============================================================
# Input and preprocess data
preprocess_data <- function(input_UMI_csv, project_name) {
    input_UMI_df <- read.csv(input_UMI_csv)
    input_UMI_matrix <- Matrix(as.matrix(input_UMI_df[, -1]))
    rownames(input_UMI_matrix) <-
        paste0(project_name, "_", input_UMI_df$Barcode)

    out <-
        #Setup the Seurat Object
        CreateSeuratObject(
            counts = t(input_UMI_matrix),
            project = project_name,
            min.cells = 3,
            min.features = 200
        ) %>%
        #Normalize the data
        NormalizeData(
            normalization.method = "LogNormalize",
            scale.factor = 10000
        ) %>%
        #Identify the highly variable features (feature selection)
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

    return(out)
}

refresh_id_in_tree <- function(integrating_tree) {
    n_hub <- nrow(integrating_tree)
    old_id <- sort(integrating_tree)
    new_id <- seq(-n_hub, n_hub) - rep(c(1, 0), c(n_hub + 1, n_hub))
    out <- matrix(
        new_id[match(integrating_tree, old_id)],
        ncol = 2, byrow = FALSE
    )
    return(out)
}
exclude_nodes_from_tree <- function(integrating_tree, excluded_node_ids) {
    excluded_rows <- c()
    for (i in seq(nrow(integrating_tree))) {
        n_excluded <- sum(integrating_tree[i, ] %in% excluded_node_ids)
        if (n_excluded == 2) {
            excluded_rows <- c(excluded_rows, i)
            excluded_node_ids <- c(excluded_node_ids, i)
        } else if (n_excluded == 1) {
            excluded_rows <- c(excluded_rows, i)
            orphan_node <- setdiff(integrating_tree[i, ], excluded_node_ids)
            integrating_tree[
                which(integrating_tree == i, arr.ind = TRUE)
            ] <- orphan_node
        }
    }
    if (length(excluded_rows) > 0) {
        integrating_tree <- integrating_tree[-excluded_rows, ]
    }
    integrating_tree <- refresh_id_in_tree(integrating_tree)
    return(integrating_tree)
}

# Integrate different SCseq samples
integrate_and_clustering <- function(
    object_list, integrating_tree, k_param
) {
    #Find integration anchors and integrate data
    integration_anchors <- FindIntegrationAnchors(
        object.list = object_list,
        anchor.features = 2000,
        scale = TRUE,
        reduction = "cca",
        l2.norm = TRUE,
        k.anchor = 5
    )

    #Run the standard workflow for visualization and clustering
    combined_object <- integration_anchors %>%
        IntegrateData(
            sample.tree = integrating_tree,
            preserve.order = TRUE,
            k.weight = 100
        ) %>%
        ScaleData() %>%
        RunPCA(npcs = 30) %>%
        FindNeighbors(reduction = "pca", dims = 1:30, k.param = k_param) %>%
        FindClusters(resolution = 0.5)

    return(combined_object)
}


# Set parameters ===============================================================
# Get input parameters from command line
input_multi_UMI_csv <- snakemake@input$multi_UMI_csv
output_integration_rds <- snakemake@output$integration_rds
integrating_order <- snakemake@params$integrating_order
k_param <- snakemake@params$k_params
n_feature_cutoff_on_sample <- snakemake@params$n_feature_cutoff_on_sample
n_barcode_cutoff_on_sample <- snakemake@params$n_barcode_cutoff_on_sample


# Implementation ===============================================================
# Create output directory
if (!dir.exists(dirname(output_integration_rds))) {
    dir.create(dirname(output_integration_rds), recursive = TRUE)
}

# Input and preprocess UMI matrix into a list
UMI_list <- sapply(
    input_multi_UMI_csv,
    function(x) {
        out <- preprocess_data(
            input_UMI_csv = x,
            project_name = file_path_sans_ext(basename(x))
        )
        return(out)
    }
)

# Remove samples with low resolution
message("Filter before integration")
is_pass <-
    sapply(
        UMI_list,
        function(UMI_object) {
            message("> ", UMI_object@project.name)
            n_feature <- UMI_object@assays$RNA@counts@Dim[1]
            n_barcode <- UMI_object@assays$RNA@counts@Dim[2]
            is_pass_feature <- (n_feature >= n_feature_cutoff_on_sample)
            is_pass_barcode <- (n_barcode >= n_barcode_cutoff_on_sample)
            message(
                ">> n_feature: ", n_feature, ": ",
                ifelse(is_pass_feature, "pass", "failed")
            )
            message(
                ">> n_barcode: ", n_barcode, ": ",
                ifelse(is_pass_barcode, "pass", "failed")
            )
            out <- is_pass_feature & is_pass_barcode
            return(out)
        }
    )
selected_UMI_list <- UMI_list[is_pass]
if (all(!is_pass)) stop("No sample pass QC.")

# Build integrating tree
integrating_tree <- matrix(integrating_order, ncol = 2, byrow = TRUE)

if (any(!is_pass)) {
    integrating_tree <-
        exclude_nodes_from_tree(
            integrating_tree,
            excluded_node_ids = -which(!is_pass)
        )
}

# Integrate UMI data
combined_object <- integrate_and_clustering(
    object_list = selected_UMI_list,
    integrating_tree = integrating_tree,
    k_param = k_param
)

# Save integrated data as rds
saveRDS(combined_object, file = output_integration_rds)
