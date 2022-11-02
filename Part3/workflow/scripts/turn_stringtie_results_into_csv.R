suppressPackageStartupMessages({
    library(dplyr)
    library(magrittr)
    library(tools)
})

# Define functions =============================================================
read_and_extract_TPM <- function(path) {
    out <-
        read.table(path, row.names = "Gene.ID",sep = "\t", header = TRUE) %>%
        select("TPM")
    return(out)
}

# Set parameters ===============================================================
# Get input parameters from snakemake
input_gene_abundance_tsvs <- snakemake@input$gene_abundance_tsv_list
output_gene_abundance_csv <- snakemake@output$gene_abundance_csv


# Implementation ===============================================================
# Input stringtie gene abundance data
gene_abundance_list <- sapply(
    simplify = FALSE,
    input_gene_abundance_tsvs,
    read_and_extract_TPM
)
names(gene_abundance_list) <-
    file_path_sans_ext(basename(input_gene_abundance_tsvs))


# Check if gene ids in files are in the same order
gene_id_list <- sapply(simplify = FALSE, gene_abundance_list, rownames)
ref_gene_id <- gene_id_list[[1]]
for (que_gene_id in gene_id_list) {
    if (length(que_gene_id) != length(ref_gene_id)) {
        stop("Files in gene_abundance_tsvs have different number of genes.")
    }
    if (!all(que_gene_id %in% ref_gene_id)) {
        stop("Files in gene_abundance_tsvs have different set of genes.")
    }
}

# Merge gene abundance of different samples
TPM_list <- sapply(
    simplify = FALSE,
    gene_abundance_list,
    magrittr::extract, ref_gene_id, "TPM"
)
gene_abundance_df <- do.call(cbind, c(list(Gene.ID = ref_gene_id), TPM_list))

write.csv(
    gene_abundance_df,
    file = output_gene_abundance_csv,
    row.names = FALSE,
    quote = FALSE
)
