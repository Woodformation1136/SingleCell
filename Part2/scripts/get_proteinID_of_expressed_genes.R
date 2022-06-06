library(Matrix)
library(tidyr)
library(dplyr)
library(magrittr)


# Define functions =============================================================
## Read UMI matrix from 10X format
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

get_express_gene_id <- function(UMI_matrix) {
    out <- rownames(UMI_matrix)[rowSums(UMI_matrix) > 0]
    return(out)
}

get_protein_id_from_gene_id <- function(gene_ids, ref_df) {
    out <- with(ref_df, protein[match(gene_ids, gene)])
    return(out)
}


# Set parameters ===============================================================
## Get input parameters from snakemake
protein2gene_csv <- snakemake@input$protein2gene_csv
input_dirs <- snakemake@input$input_dirs
expressed_protein_id_txt <- snakemake@output$expressed_protein_id_txt

# protein2gene_csv <- "results/protein2gene.csv"
# output_txt <- "results/Expressed_proteinID_list/Species4_Sample5.txt"
# input_dirs <- c(
#     "results/SCseq_TenX_Ptr",
#     "results/SCseq_TenX_PalChen2021",
#     "results/SCseq_TenX_Lch",
#     "results/SCseq_MARSseq_Egr",
#     "results/SCseq_MARSseq_Tar"
# )

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


# Implementation ===============================================================
# Input data
input_UMI_matrix <- input_matrix_dirs %>%
    sapply(read_UMI_matrix, simplify = FALSE) %>%
    set_names(input_samples)
protein2gene_df <- read.csv(protein2gene_csv)
stopifnot(length(unique(protein2gene_df$protein)) == nrow(protein2gene_df))
stopifnot(length(unique(protein2gene_df$gene)) == nrow(protein2gene_df))

# Get expressed gene ID
expressed_gene_id <- sapply(
    input_UMI_matrix,
    get_express_gene_id
)

# Convert gene ID into protein ID
expressed_protein_id <- sapply(
    expressed_gene_id,
    get_protein_id_from_gene_id,
    ref_df = protein2gene_df
)

# Output expressed protein ID
writeLines(
    unlist(expressed_protein_id),
    expressed_protein_id_txt
)