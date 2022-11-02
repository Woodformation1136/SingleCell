suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
})

# Define functions =============================================================


# Set parameters ===============================================================
# Get input parameters from snakemake
input_SC_geneUMI_csv <- snakemake@input$geneUMI_csv
input_SC_plotting_csv <- snakemake@input$plotting_csv
input_LCM_gene_abundance_csv <- snakemake@input$gene_abundance_csv
output_plotting_csv <- snakemake@output$plotting_csv
cor_method <- snakemake@params$cor_method


# Implementation ===============================================================
# Input data
SC_geneUMI_df <- read.csv(input_SC_geneUMI_csv, row.names = "Barcode")
SC_plotting_df <- read.csv(input_SC_plotting_csv, row.names = "Barcode")
LCM_gene_abundance_df <- read.csv(
    input_LCM_gene_abundance_csv,
    row.names = "Gene.ID"
)

# Calculate mean gene abundance of LCM data
LCM_gene_abundance_df$Mean_gene_abundances <-
    apply(LCM_gene_abundance_df, 1, mean)

# Calculate correlation coefficients between LCM samples and SC barcodes
ref_gene_id <- colnames(SC_geneUMI_df)
cor_matrix <- cor(
    t(SC_geneUMI_df),
    LCM_gene_abundance_df[ref_gene_id, ],
    method = cor_method
)
colnames(cor_matrix) <-
    colnames(cor_matrix) %>%
    sub(pattern = "_gene_abundances", replacement = "") %>%
    paste0("Correlation_", .)

# Merge correlation coefficients into plotting data frame
stopifnot(all(rownames(SC_plotting_df) == rownames(cor_matrix)))
merged_plotting_df <- cbind(
    Barcode = rownames(SC_plotting_df), SC_plotting_df, cor_matrix
)

write.csv(
    merged_plotting_df, file = output_plotting_csv,
    row.names = FALSE, quote = FALSE
)
