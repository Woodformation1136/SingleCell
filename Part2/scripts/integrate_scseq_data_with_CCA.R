library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)


# Define functions =============================================================
# Input and preprocess 10X data
preprocess_10X_data <- function(data_dir, project_name) {
    out <- Read10X(data.dir = data_dir) %>%
    #Setup the Seurat Object
    CreateSeuratObject(
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

# Integrate different SCseq samples
integrate_and_clustering <- function(object_list, integrating_tree, k_param) {
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
            preserve.order = TRUE
        ) %>%
        ScaleData() %>%
        RunPCA(npcs = 30) %>%
        RunUMAP(reduction = "pca", dims = 1:30) %>%
        FindNeighbors(reduction = "pca", dims = 1:30, k.param = k_param) %>%
        FindClusters(resolution = 0.5)

    return(combined_object)
}


# Set parameters ===============================================================
# Get input parameters from command line
input_expression_dir <- snakemake@input$input_expression_dir
output_path <- snakemake@output$output_path
k_param <- snakemake@params$k_param
integration_order <- snakemake@params$integration_order
# input_expression_dir = "results/Ortholog_UMI/Species4_Sample5"
# output_path = "results/Integrated_rds/temp.rds"
# integration_order = c(
#     "10X_Ptr", "10X_Ptr_Chen2021", "10X_Lch", "MARSseq_Egr", "MARSseq_Tar"
# )


# Implementation ===============================================================
# Create output directory
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input 10X data
## Input and preprocess transcriptome data
exp_object <- sapply(
    integration_order,
    function(x) {
        exp_path <- paste0(input_expression_dir, "/", x)
        out <- preprocess_10X_data(
            data_dir = exp_path,
            project_name = x
        )
        return(out)
    }
)

# Integrate transcriptome data
integrating_index <- (-match(integration_order, names(exp_object)))
if (length(integrating_index) == 2) {
    integrating_tree <- matrix(integrating_index, byrow = TRUE, ncol = 2)
} else {
    n_node <- length(integrating_index) - 2
    integrating_tree <- cbind(
        c(integrating_index[1], seq(n_node)),
        integrating_index[-1]
    )
}
combined_object <- integrate_and_clustering(
    object_list = exp_object,
    integrating_tree = integrating_tree,
    k_param = k_param
)

# Save integrated data as rds
saveRDS(combined_object, file = output_path)
