library(magrittr)
library(Matrix)
library(dplyr)


# Define functions =============================================================
# Calculate the ortholog transcript abundance
calculate_ortholog_UMI_matrix <-
  function(UMI_matrix,
           ortholog_table,
           selected_ortholog_cluster) {
    matrix_geneID <- rownames(UMI_matrix)
    geneID_list_of_each_cluster <-
      sapply(
        selected_ortholog_cluster,
        function(cc) {
          with(ortholog_table,
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
ortholog_long_csv <- snakemake@input$ortholog_long_convertedID_csv
input_geneUMI_csv <- snakemake@input$geneUMI_csv
output_orthologUMI_csv <- snakemake@output$orthologUMI_csv


# Implementation ===============================================================
# Create output directory
if (!dir.exists(dirname(output_orthologUMI_csv))) {
  dir.create(dirname(output_orthologUMI_csv), recursive = TRUE)
}

# Input data
ortholog_long <- read.csv(ortholog_long_csv)

input_geneUMI_df <- read.csv(input_geneUMI_csv)
input_geneUMI_matrix <- Matrix(as.matrix(input_geneUMI_df[, -1]))
rownames(input_geneUMI_matrix) <- input_geneUMI_df$Barcode

# Select ortholog groups containing gene ID within gene UMI count matrix
selected_ortholog_cluster <-
  ortholog_long %>%
    filter(gene_name %in% colnames(input_geneUMI_matrix)) %>%
    extract(, "cluster_id") %>%
    unique()

# Calculate ortholog UMI matrix
output_orthologUMI_matrix <-
  calculate_ortholog_UMI_matrix(
    UMI_matrix = t(input_geneUMI_matrix),
    ortholog_table = ortholog_long,
    selected_ortholog_cluster = selected_ortholog_cluster
  )

# Create and output ortholog UMI count table
output_orthologUMI_df <- cbind(
    Barcode = colnames(output_orthologUMI_matrix),
    data.frame(as.matrix(t(output_orthologUMI_matrix)))
)

write.csv(
  output_orthologUMI_df,
  file = output_orthologUMI_csv,
  row.names = FALSE,
  quote = FALSE
)
