library(magrittr)
library(getopt)
library(Matrix)
library(dplyr)


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

# Construct the ortholog transcript abundance matrix for SC data
# We use the preprocessed data above to construct the ortholog UMI matrix based
# on the overlap ortholog clusters among certain species, and output as the
# format of CellRanger output (we called it 10X format below).

## Calculate the ortholog transcript abundance
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

## Output ortholog UMI matrix as 10X matrix format
output_as_10X_format <- function(matrix, output_dir, barcode_prefix = "Ptr") {
  gz_temp <- gzfile(paste0(output_dir, "features.tsv.gz"), "w")
  feature_names <- rownames(matrix)
  write.table(
    data.frame(feature_names, feature_names, "Gene Expression"),
    gz_temp,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  close(gz_temp)

  gz_temp <- gzfile(paste0(output_dir, "barcodes.tsv.gz"), "w")
  write.table(
    data.frame(paste0(barcode_prefix, colnames(matrix))),
    gz_temp,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  close(gz_temp)

  writeMM(matrix, paste0(output_dir, "matrix.mtx"))
  system(
    paste0("gzip --force ", paste0(output_dir, "matrix.mtx"))
  )

  return("Finish")
}


# Set parameters ===============================================================
# Get input parameters from command line
# spec <- matrix(
#   c(
#     "input_dirs", "m", 2, "character", "Folders from cellranger reanalysis",
#     "ortholog_long_csv", "i", 2, "character", "Long ortholog table",
#     "output_dir", "o", 2, "character", "Output folder"
#   ),
#   byrow = TRUE, ncol = 5
# )
# opt <- getopt(spec = spec)
# input_dirs <- strsplit(opt$input_dirs, split = ",")[[1]]
# ortholog_long_csv <- opt$ortholog_long_csv
# output_dir <- opt$output_dir

input_dirs <- snakemake@input$input_dirs
ortholog_long_csv <- snakemake@input$ortholog_long_csv
output_dir <- snakemake@output$output_dir

input_samples <- sapply(
  input_dirs,
  function(x) {
    suffix <- sub(pattern = "results/SCseq_", replacement = "", x)
    return(suffix)
  }
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
## Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
for (dir in input_samples) {
  if (!dir.exists(paste(output_dir, dir, sep = "/"))) {
    dir.create(paste(output_dir, dir, sep = "/"))
  }
}

## Input data
input_UMI_matrix <- sapply(input_matrix_dirs, read_UMI_matrix, simplify = FALSE)
names(input_UMI_matrix) <- input_samples
ortholog_long <- read.csv(ortholog_long_csv)


## If gene ids are consistent among samples, directly output UMI matrix;
## otherwise, calculate ortholog UMI matrix before output.
all_gene_id  <- c()
for (i in seq_along(input_UMI_matrix)) {
  all_gene_id <- union(all_gene_id, rownames(input_UMI_matrix[[i]]))
}

if (nrow(input_UMI_matrix[[1]]) == length(all_gene_id)) {
  for (i in seq_along(input_UMI_matrix)) {
    ## Output UMI matrix
      suffix <- names(input_UMI_matrix)[i]
      output_as_10X_format(
        input_UMI_matrix[[i]][all_gene_id, ],
        output_dir = paste0(output_dir, "/", suffix, "/"),
        barcode_prefix = paste0(suffix, "_")
      )
  }
} else {
  ## Select ortholog groups (overlap_clusters)
  for (i in seq_along(input_UMI_matrix)) {
    if (i == 1) {
      overlap_clusters <-
        ortholog_long %>%
        filter(gene_name %in% rownames(input_UMI_matrix[[i]])) %>%
        extract(, "cluster_id")
    } else {
      overlap_clusters <-
        intersect(
          overlap_clusters,
          ortholog_long %>%
            filter(gene_name %in% rownames(input_UMI_matrix[[i]])) %>%
            extract(, "cluster_id")
        )
    }
  }

  ## Calculate and output ortholog UMI matrix
  for (i in seq_along(input_UMI_matrix)) {
    ortho_UMI_matrix <-
      calculate_ortholog_UMI_matrix(
        UMI_matrix = input_UMI_matrix[[i]],
        ortholog_table = ortholog_long,
        selected_ortholog_cluster = overlap_clusters
      )

    suffix <- names(input_UMI_matrix)[i]
    output_as_10X_format(
      ortho_UMI_matrix,
      output_dir = paste0(output_dir, "/", suffix, "/"),
      barcode_prefix = paste0(suffix, "_")
    )
  }
}
