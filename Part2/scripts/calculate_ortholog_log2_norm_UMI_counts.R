library(magrittr)
library(dplyr)
library(RColorBrewer)
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

# Read UMI matrix from 10X format
read_UMI_matrix <- function(matrix_dir) {
    barcode_path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features_path <- paste0(matrix_dir, "features.tsv.gz")
    matrix_path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix_path)
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
    colnames(mat) <- barcode_names$V1
    rownames(mat) <- feature_names$V1
    return(mat)
}

# Calculate normalized UMI in log2 scale
calcualte_log_normalized_UMI <- function(UMI_matrix) {
    norm_factor <- 1000 / colSums(UMI_matrix)
    norm_UMI_matrix <- apply(UMI_matrix, 1, multiply_by, norm_factor) %>% t()
    log2_norm_UMI_matrix <- log2(norm_UMI_matrix + 1)
    return(log2_norm_UMI_matrix)
}

# Calculate the ortholog transcript abundance
calculate_ortholog_UMI_matrix <- function(
    UMI_matrix,
    ortholog_table,
    selected_ortholog_cluster
) {
    matrix_geneID <- rownames(UMI_matrix)
    geneID_list_of_each_cluster <-
        sapply(
            selected_ortholog_cluster,
            function(cc) {
                with(
                    ortholog_table,
                    intersect(matrix_geneID, gene_name[cluster_id == cc])
                )
            }
        )
    cluster_indicator_matrix <-
        Matrix(0,
            nrow = length(selected_ortholog_cluster),
            ncol = length(matrix_geneID),
            dimnames = list(
                paste0("Cluster_", selected_ortholog_cluster),
                matrix_geneID
            )
        )
    for (ci in seq_along(geneID_list_of_each_cluster)) {
        selected_genes <- geneID_list_of_each_cluster[[ci]]
        cluster_indicator_matrix[ci, selected_genes] <- 1
    }
    ortho_UMI_matrix <- cluster_indicator_matrix %*% UMI_matrix
    return(ortho_UMI_matrix)
}


# Set parameters ===============================================================
# Get input parameters from snakemake
input_dirs <- snakemake@input$input_dirs
ortholog_long_csv <- snakemake@input$ortholog_long_csv
output_dir <- snakemake@output$output_dir

# input_dirs <- c(
#     "results/SCseq_TenX_Ptr",
#     "results/SCseq_TenX_PalChen2021",
#     "results/SCseq_TenX_Lch",
#     "results/SCseq_MARSseq_Egr",
#     "results/SCseq_MARSseq_Tar"
# )
# ortholog_long_csv <- "results/Ortholog_table/all_group_long.csv"

input_samples <- sapply(
    input_dirs,
    sub,
    pattern = "results/SCseq_",
    replacement = ""
)

input_matrix_dirs <- sapply(
    input_samples,
    function(x) {
        out <- paste0(
            "results/SCseq_", x, "/cellranger_reanalysis_", x,
            "/outs/filtered_feature_bc_matrix/"
        )
        return(out)
    }
)

# Implement program ============================================================
# Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input data
input_UMI_matrix <- input_matrix_dirs %>%
    sapply(read_UMI_matrix, simplify = FALSE) %>%
    set_names(input_samples)
ortho_table <- read.csv(ortholog_long_csv)

# Output log2-transformed normalized ortholog UMI matrix for each sample
for (sample_name in input_samples) {
    # Select ortholog groups (overlap_clusters)
    selected_gene_name <-
        input_UMI_matrix[[sample_name]] %>%
        rownames() %>%
        unlist()
    selected_cluster_id <-
        ortho_table %>%
        filter(gene_name %in% selected_gene_name) %>%
        extract(, "cluster_id") %>%
        unique()

    # Calculate normalized ortholog UMI in log2 scale
    log2_norm_ortho_UMI_matrix <-
        input_UMI_matrix[[sample_name]] %>%
        calcualte_log_normalized_UMI() %>%
        calculate_ortholog_UMI_matrix(
            ortholog_table = ortho_table,
            selected_ortholog_cluster = selected_cluster_id
        )

    # Output matrix as rds
    saveRDS(
        object = log2_norm_ortho_UMI_matrix,
        file = paste0(output_dir, "/", sample_name, ".rds")
    )
}
