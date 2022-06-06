library(magrittr)
library(Matrix)
library(dplyr)

setwd("C:/Users/ftroo/Desktop/SC_LCM/00_reorganize_the_data")

# Load and process the gene transcript abundance data ========================
# Load in the following data:
# 1. LCM analysis results from StringTie and DESeq2
# 2. SC analysis results from CellRanger and MARS-seq2.0
# 3. ortholog clustering results from orthoMCL
# Then reorganize and process the data into matrix or data frame.
dir.create("output_RDS")
dir.create("output_10X_matrix")

# Load and process LCM data ----------------------------------------------------
## Define the loading function
load_LCM_data_into_matrix <- function(TPM_path, label_table) {
  # Input the label table
  label_df <- read.csv(label_table)
  # Input and combine the TPM of each sample
  for (i in seq(nrow(label_df))) {
    TPM_column <-
      read.table(paste0(TPM_path, "/", label_df$Group[i], ".tsv"),
                 row.names = "Gene.ID", sep = "\t", header = TRUE) %>%
      extract(, "TPM", drop = FALSE)
    if (i == 1) {
      TPM_matrix <- TPM_column
      gene_ID <- rownames(TPM_column)
    } else {
      TPM_matrix <- cbind(TPM_matrix, TPM_column[gene_ID, ])
    }
  }
  colnames(TPM_matrix) <- label_df$Condition
  TPM_matrix <- Matrix(as.matrix(TPM_matrix))
  TPM_matrix
}

## Define the function to calculate the mean TPM for each cell type
calculate_mean_TPM <- function(TPM_matrix, cell_type) {
  mean_TPM_matrix <-
    sapply(cell_type,
           function(cc) {
             selected_columns <- grepl(cc, colnames(TPM_matrix))
             mean_TPM_column <- apply(TPM_matrix[, selected_columns], 1, mean)
             mean_TPM_column
           }
    )
  mean_TPM_matrix <- Matrix(mean_TPM_matrix)
  mean_TPM_matrix
}

## Define the function to extract DEGs for each cell type
## (Adjusted P values have already been filtered in previous analysis pipeline)
extract_LCM_DEG_list <- function(DEG_table_prefix, cell_type, log2FC) {
  LCM_DEG_list <- list()
  for (i in seq_along(cell_type)) {
    LCM_DEG_table <- read.csv(paste0(DEG_table_prefix, cell_type[i], ".csv"))
    log2FC_column_index <- grepl("log2FoldChange", colnames(LCM_DEG_table))
    is_up_reg <-
      LCM_DEG_table[, log2FC_column_index, drop = F] %>%
      apply(1, function(X) all(X >= log2FC))
    is_down_reg <-
      LCM_DEG_table[, log2FC_column_index, drop = F] %>%
      apply(1, function(X) all(X <= -log2FC))
    LCM_DEG_list[[i]] <- list(
      up_reg_genes = LCM_DEG_table$X[is_up_reg],
      down_reg_genes = LCM_DEG_table$X[is_down_reg]
    )
  }
  names(LCM_DEG_list) <- cell_type
  LCM_DEG_list
}

## Load and process Ptr LCM data
Ptr_LCM_TPM_matrix <-
  load_LCM_data_into_matrix("Ptr_LCM",
                            "Ptr_LCM/20210127_Inputfilelist.csv")
Ptr_LCM_mean_TPM_matrix <-
  calculate_mean_TPM(Ptr_LCM_TPM_matrix,
                     c("vessel", "ray", "fiber"))
Ptr_LCM_DEG_list <-
  extract_LCM_DEG_list("Ptr_LCM/outputDESeq2_DEGtable_Populus_trichocarpa_",
                       c("vessel", "ray", "fiber"),
                       log2FC = 1)

saveRDS(Ptr_LCM_TPM_matrix,
        "output_RDS/RDS_Ptr_LCM_TPM_matrix.rds")
saveRDS(Ptr_LCM_mean_TPM_matrix,
        "output_RDS/RDS_Ptr_LCM_mean_TPM_matrix.rds")
saveRDS(Ptr_LCM_DEG_list,
        "output_RDS/RDS_Ptr_LCM_DEG_list.rds")

# Load and process SC data -----------------------------------------------------
## Define the loading function for CellRanger
load_SC_data_into_matrix_CellRanger <- function(matrix_path) {
  mat <- readMM(file = paste0(matrix_path, "/matrix.mtx.gz"))
  feature.names <- read.delim(paste0(matrix_path, "/features.tsv.gz"),
                              header = FALSE,
                              stringsAsFactors = FALSE)
  barcode.names <- read.delim(paste0(matrix_path,"/barcodes.tsv.gz"),
                              header = FALSE,
                              stringsAsFactors = FALSE)
  colnames(mat) <- barcode.names$V1
  rownames(mat) <- feature.names$V1
  SCseq_UMI_Matrix <- Matrix(mat, sparse = TRUE)
  SCseq_UMI_Matrix
}

## Define the loading function for MARS-seq2.0
load_SC_data_into_matrix_MARS <- function(CDS_list, UMI_table_path) {
  selected_CDS <- scan(CDS_list, what = "")
  file_list <- list.files(UMI_table_path)
  UMI_matrix <- 
    read.table(paste0(UMI_table_path, "/", file_list[1]),
               header = TRUE, sep = "\t") %>%
    extract(selected_CDS, ) %>%
    as.matrix() %>%
    Matrix(sparse = TRUE)
  cat("\r", 1, "/",length(file_list))
    
  for (f in 2:length(file_list)) {
    foo <- read.table(paste0(UMI_table_path, "/", file_list[f]),
                      header = TRUE, sep = "\t") %>%
      extract(selected_CDS, ) %>%
      as.matrix() %>%
      Matrix(sparse = TRUE)
    UMI_matrix <- cbind(UMI_matrix, foo)
    cat("\r", f, "/",length(file_list))
  }
  UMI_matrix
}

## Define the function to calculate the mean counts for each cell cluster
calculate_mean_counts <- function(UMI_matrix,
                                  cluster_table) {
  n_clusters <- max(cluster_table$Cluster)
  mean_count_matrix <-
    sapply(1:n_clusters, function(cl) {
      selected_barcodes <- with(cluster_table, Barcode[Cluster == cl])
      mean_count_column <-
        UMI_matrix %>%
        extract(, selected_barcodes, drop = FALSE) %>%
        apply(1, mean)
      mean_count_column
    })
  colnames(mean_count_matrix) <-
    paste0("CellCluster_", 1:n_clusters)
  mean_count_matrix <- Matrix(mean_count_matrix, sparse = TRUE)
  mean_count_matrix
}

## Define the function to extract DEGs for each cluster
extract_SC_DEG_list <- function(SC_DEG_table, log2FC, FDR) {
  n_clusters <- (ncol(SC_DEG_table) - 2) / 3
  SC_DEG_list <- 
    sapply(seq(n_clusters),
           simplify = F,
           function(cc) {
             is_up_reg <- (SC_DEG_table[, 3 * cc + 1] >= log2FC &
                           SC_DEG_table[, 3 * cc + 2] < FDR)
             is_down_reg <- (SC_DEG_table[, 3 * cc + 1] <= -log2FC &
                             SC_DEG_table[, 3 * cc + 2] < FDR)
             DEG <- list(
               up_reg_genes = SC_DEG_table$Feature.ID[is_up_reg],
               down_reg_genes = SC_DEG_table$Feature.ID[is_down_reg]
             )
             DEG
           }
  )
  names(SC_DEG_list) <- paste0("Cluster_", seq(n_clusters))
  SC_DEG_list
}

## Load and process Ptr SC data
Ptr_SC_UMI_matrix <-
  load_SC_data_into_matrix_CellRanger(
    "Ptr_SC_CellRanger/filtered_feature_bc_matrix"
  )
Ptr_SC_norm_factor <- 1000 / colSums(Ptr_SC_UMI_matrix)
Ptr_SC_norm_UMI_matrix <-
  Ptr_SC_UMI_matrix %>%
  apply(1, multiply_by, Ptr_SC_norm_factor) %>%
  Matrix() %>% t()
Ptr_SC_cell_cluster_10 <-
  read.csv("Ptr_SC_CellRanger/clustering/kmeans_10_clusters/clusters.csv")
Ptr_SC_mean_norm_UMI_matrix <-
  calculate_mean_counts(Ptr_SC_norm_UMI_matrix,
                        Ptr_SC_cell_cluster_10)
Ptr_SC_projection_umap <-
  read.csv("Ptr_SC_CellRanger/umap/2_components/projection.csv")
Ptr_SCseq_DEG_table <-
  read.csv(
    "Ptr_SC_CellRanger/diffexp/kmeans_10_clusters/differential_expression.csv"
  )
Ptr_SCseq_DEG_list <-
  extract_SC_DEG_list(Ptr_SCseq_DEG_table, log2FC = 1, FDR = 0.05)

saveRDS(Ptr_SC_UMI_matrix,
        "output_RDS/RDS_Ptr_SC_UMI_matrix.rds")
saveRDS(Ptr_SC_norm_factor,
        "output_RDS/RDS_Ptr_SC_norm_factor.rds")
saveRDS(Ptr_SC_norm_UMI_matrix,
        "output_RDS/RDS_Ptr_SC_norm_UMI_matrix.rds")
saveRDS(Ptr_SC_cell_cluster_10,
        "output_RDS/RDS_Ptr_SC_cell_cluster_10.rds")
saveRDS(Ptr_SC_mean_norm_UMI_matrix,
        "output_RDS/RDS_Ptr_SC_mean_norm_UMI_matrix.rds")
saveRDS(Ptr_SC_projection_umap,
        "output_RDS/RDS_Ptr_SC_projection_umap.rds")
saveRDS(Ptr_SCseq_DEG_table,
        "output_RDS/RDS_Ptr_SCseq_DEG_table.rds")
saveRDS(Ptr_SCseq_DEG_list,
        "output_RDS/RDS_Ptr_SCseq_DEG_list.rds")

## Load and process Lch SC data
Lch_SC_UMI_matrix <-
  load_SC_data_into_matrix_CellRanger(
    "Lch_SC_CellRanger/filtered_feature_bc_matrix"
  )
Lch_SC_norm_factor <- 1000 / colSums(Lch_SC_UMI_matrix)
Lch_SC_cell_cluster_10 <-
  read.csv("Lch_SC_CellRanger/clustering/kmeans_10_clusters/clusters.csv")

saveRDS(Lch_SC_UMI_matrix,
        "output_RDS/RDS_Lch_SC_UMI_matrix.rds")
saveRDS(Lch_SC_norm_factor,
        "output_RDS/RDS_Lch_SC_norm_factor.rds")
saveRDS(Lch_SC_cell_cluster_10,
        "output_RDS/RDS_Lch_SC_cell_cluster_10.rds")

## Load and process Egr SC data
Egr_SC_UMI_matrix <-
  load_SC_data_into_matrix_MARS(
    CDS_list = "MARS-seq_UMItable/Eg_CDS.txt",
    UMI_table_path = "MARS-seq_UMItable/MARS-seq_UMItable_Eg")
selected_barcodes <- (colSums(Egr_SC_UMI_matrix) >= 100)
Egr_SC_UMI_matrix <- Egr_SC_UMI_matrix[, selected_barcodes]
Egr_SC_norm_factor <- 1000 / colSums(Egr_SC_UMI_matrix)
Egr_SC_cell_cluster_10 <-
  read.csv("Egr_SC_CellRanger/clustering/kmeans_10_clusters/clusters.csv")

saveRDS(Egr_SC_UMI_matrix,
        "output_RDS/RDS_Egr_SC_UMI_matrix.rds")
saveRDS(Egr_SC_norm_factor,
        "output_RDS/RDS_Egr_SC_norm_factor.rds")
saveRDS(Egr_SC_cell_cluster_10,
        "output_RDS/RDS_Egr_SC_cell_cluster_10.rds")

## Load and process Egr SC data
Tar_SC_UMI_matrix <-
  load_SC_data_into_matrix_MARS(
    CDS_list = "MARS-seq_UMItable/Ta_CDS.txt",
    UMI_table_path = "MARS-seq_UMItable/MARS-seq_UMItable_Ta")
rownames(Tar_SC_UMI_matrix) <-
  paste0('evm.TU.',rownames(Tar_SC_UMI_matrix))
selected_barcodes <- (colSums(Tar_SC_UMI_matrix) >= 100)
Tar_SC_UMI_matrix <- Tar_SC_UMI_matrix[, selected_barcodes]
Tar_SC_norm_factor <- 1000 / colSums(Tar_SC_UMI_matrix)
Tar_SC_cell_cluster_18 <-
  read.csv("Tar_SC_CellRanger/clustering/kmeans_18_clusters/clusters.csv")
Tar_SC_cell_cluster_18$Barcode <-
  sub("Kun", "Ta", Tar_SC_cell_cluster_18$Barcode)

saveRDS(Tar_SC_UMI_matrix,
        "output_RDS/RDS_Tar_SC_UMI_matrix.rds")
saveRDS(Tar_SC_norm_factor,
        "output_RDS/RDS_Tar_SC_norm_factor.rds")
saveRDS(Tar_SC_cell_cluster_18,
        "output_RDS/RDS_Tar_SC_cell_cluster_18.rds")

# Load and process the ortholog data -------------------------------------------
## Load and covert the ortholog table from wide data to long data
cat(file = "Geneclustering_OrthoMCL/all_group_long.csv",
    paste("clusterID","speciesName","geneName\n", sep = ","))

ortholog_table_wide <- file("Geneclustering_OrthoMCL/all.group", "r")
reach_EOF <- FALSE
iterator_index <- 1
while (!reach_EOF) {
  ortholog_row <- readLines(ortholog_table_wide, n = 1)
  reach_EOF <- (length(ortholog_row) == 0)
  if (reach_EOF) break()
  
  cut_row <- ortholog_row %>%
    gsub(pattern = '"', replacement = "") %>%
    strsplit(split = "\t") %>%
    extract2(1)
  clusterID <- cut_row[1] %>%
    strsplit(split = "[(]") %>%
    extract2(1) %>%
    extract(1)
  species_gene <- cut_row[-1] %>%
    strsplit(split = "[(]") %>%
    sapply(extract, 1)
  stopifnot(all(regexpr("_", species_gene) == 4)) # TRUE
  speciesName <- substr(species_gene, 1, 3)
  geneName <- substring(species_gene, 5)
  
  outDataFrame <- data.frame(clusterID, speciesName, geneName)
  
  write.table(outDataFrame,
              file = "Geneclustering_OrthoMCL/all_group_long.csv",
              row.names = FALSE, col.names = FALSE,
              append = TRUE, quote = FALSE, sep = ",")

  cat("\r", iterator_index, "/ 133114")
  iterator_index %<>% add(1)
}
close(ortholog_table_wide)

## Load the ortholog table (long data format)
ortho_table <- read.csv("Geneclustering_OrthoMCL/all_group_long.csv")
dim(ortho_table) # 478525 3

## Convert the protein ID into gene ID for Ptr
### Define the function for ID extracting
extract_prefix_before_the_point <-
  function(raw_ID,
           which_point = 2) {
    cut_position <-
      raw_ID %>%
      gregexpr(pattern = "[.]") %>%
      extract2(1) %>% 
      extract(which_point)
    extracted_ID <- substr(raw_ID, 1, cut_position - 1)
    extracted_ID
  }
### Replace Ptr peptide ID with gene ID
Ptr_protein_fasta <- 
  scan("Ptrichocarpa_533_v4.1.protein_primaryTranscriptOnly.fa",
       what = "", sep = "\n")
fasta_info <- grep(">", Ptr_protein_fasta, value = TRUE)
Ptr_protein_ID <-
  strsplit(fasta_info, " ") %>%
  sapply(extract, 1) %>%
  sub(">", "", .)
Ptr_gene_ID <-
  strsplit(fasta_info, " ") %>%
  sapply(extract, 5) %>%
  sub("ID=", "", .) %>%
  sapply(extract_prefix_before_the_point, 2)
stopifnot(
  all(filter(ortho_table, speciesName == "PoT")$geneName %in% Ptr_protein_ID)
)
stopifnot(
  all(rownames(Ptr_SC_norm_UMI_matrix) %in% paste0(Ptr_gene_ID, ".v4.1"))
)
Ptr_row_indices = which(ortho_table$speciesName == "PoT")
ortho_table[Ptr_row_indices, "geneName"] <-
  paste0(Ptr_gene_ID[match(ortho_table[Ptr_row_indices, "geneName"],
                           Ptr_protein_ID)], ".v4.1")
stopifnot(
  all(rownames(Ptr_SC_norm_UMI_matrix) %in% ortho_table$geneName)
)

### Replace Egr peptide ID with gene ID
Egr_row_indices = which(ortho_table$speciesName == "EuG")
ortho_table[Egr_row_indices, "geneName"] <-
  sub("1.p", "v2.0", ortho_table[Egr_row_indices, "geneName"])

### Replace Tar peptide ID with gene ID
Tar_gtf <-
  read.table("Trochodendron_aralioides_chromosomes_pasa2.longest.filter.gtf",
             sep = "\t", header = FALSE)
ID_info <- unique(Tar_gtf$V9)
Tar_transcript_id <-
  strsplit(ID_info, ";") %>%
  sapply(extract, 1) %>%
  sub("transcript_id ", "", .)
Tar_gene_id <-
  strsplit(ID_info, ";") %>%
  sapply(extract, 2) %>%
  sub(" gene_id ", "", .)
stopifnot(
  all(filter(ortho_table, speciesName == "TrA")$geneName %in% Tar_transcript_id)
)
stopifnot(
  all(rownames(Tar_SC_UMI_matrix) %in% Tar_gene_id)
)
Tar_row_indices <- which(ortho_table$speciesName == "TrA")
ortho_table[Tar_row_indices, "geneName"] <-
  Tar_gene_id[match(ortho_table[Tar_row_indices, "geneName"],
                    Tar_transcript_id)]
stopifnot(
  all(rownames(Tar_SC_UMI_matrix) %in% ortho_table$geneName)
)

# saveRDS(ortho_table, "RDS_ortho_cluster_table.rds")

Ptr_ortho_table = filter(ortho_table, speciesName == "PoT")
Egr_ortho_table = filter(ortho_table, speciesName == "EuG")
Tar_ortho_table = filter(ortho_table, speciesName == "TrA")
Lch_ortho_table = filter(ortho_table, speciesName == "LiC")

saveRDS(Ptr_ortho_table,
        "output_RDS/RDS_Ptr_ortho_table.rds")
saveRDS(Egr_ortho_table,
        "output_RDS/RDS_Egr_ortho_table.rds")
saveRDS(Tar_ortho_table,
        "output_RDS/RDS_Tar_ortho_table.rds")
saveRDS(Lch_ortho_table,
        "output_RDS/RDS_Lch_ortho_table.rds")

# Construct the ortholog transcript abundance matrix for SC data ===============
# We use the preprocessed data above to construct the ortholog UMI matrix based
# on the overlap ortholog clusters among certain species, and output as the
# format of CellRanger output (we called it 10X format below).

## Define the function to calculate the ortholog transcript abundance
calculate_ortholog_UMI_matrix <-
  function(UMI_matrix,
           ortholog_table,
           selected_ortholog_cluster){
    geneID_list_of_each_cluster <-
      sapply(selected_ortholog_cluster,
             function(cc) with(ortholog_table, geneName[clusterID == cc]))
    matrix_geneID <- rownames(UMI_matrix)
    cluster_indicator_matrix <-
      Matrix(0,
             nrow = length(selected_ortholog_cluster),
             ncol = length(matrix_geneID),
             dimnames = list(
               paste0("Cluster_", selected_ortholog_cluster),
               matrix_geneID
             )
      )
    for (i in seq_along(geneID_list_of_each_cluster)) {
      selected_genes <- geneID_list_of_each_cluster[[i]]
      cluster_indicator_matrix[i, selected_genes] <- 1
    }
    ortho_UMI_matrix <- cluster_indicator_matrix %*% UMI_matrix
    ortho_UMI_matrix
}

## Define the function to output the ortholog UMI matrix as 10X matrix format
output_as_10X_format <- function(matrix, prefix) {
  gz_temp <- gzfile(paste0(prefix, "features.tsv.gz"), "w")
  foo <- rownames(matrix)
  write.table(
    data.frame(foo, foo, "Gene Expression"),
    gz_temp,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  close(gz_temp)

  gz_temp <- gzfile(paste0(prefix, "barcodes.tsv.gz"), "w")
  write.table(
    data.frame(paste0(prefix, colnames(matrix))),
    gz_temp,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  close(gz_temp)

  writeMM(matrix, paste0(prefix, "matrix.mtx"))
  system(paste0("gzip ", paste0(prefix, "matrix.mtx")))
  
  "Finish"
}

## Construct the ortholog UMI matrix of 4 species (Ptr, Egr, Tar and Lch)
overlap_clusters <- intersect(intersect(unique(Ptr_ortho_table$clusterID),
                                        unique(Egr_ortho_table$clusterID)),
                              intersect(unique(Tar_ortho_table$clusterID),
                                        unique(Lch_ortho_table$clusterID)))
Ptr_SC_4S_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Ptr_SC_UMI_matrix,
    Ptr_ortho_table,
    overlap_clusters
  )
Egr_SC_4S_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Egr_SC_UMI_matrix,
    Egr_ortho_table,
    overlap_clusters
  )
Tar_SC_4S_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Tar_SC_UMI_matrix,
    Tar_ortho_table,
    overlap_clusters
  )
Lch_SC_4S_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Lch_SC_UMI_matrix,
    Lch_ortho_table,
    overlap_clusters
  )

output_as_10X_format(Ptr_SC_4S_ortho_UMI_matrix,
                     "output_10X_matrix/RDS_Ptr_SC_4S_ortho_UMI_matrix.rds")
output_as_10X_format(Egr_SC_4S_ortho_UMI_matrix,
                     "output_10X_matrix/RDS_Egr_SC_4S_ortho_UMI_matrix.rds")
output_as_10X_format(Tar_SC_4S_ortho_UMI_matrix,
                     "output_10X_matrix/RDS_Tar_SC_4S_ortho_UMI_matrix.rds")
output_as_10X_format(Lch_SC_4S_ortho_UMI_matrix,
                     "output_10X_matrix/RDS_Lch_SC_4S_ortho_UMI_matrix.rds")


## Construct the ortholog UMI matrix of Ptr and Egr
overlap_clusters <- intersect(unique(Ptr_ortho_table$clusterID),
                              unique(Egr_ortho_table$clusterID))
Ptr_SC_PtrEgr_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Ptr_SC_UMI_matrix,
    Ptr_ortho_table,
    overlap_clusters
  )
Egr_SC_PtrEgr_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Egr_SC_UMI_matrix,
    Egr_ortho_table,
    overlap_clusters
  )

Ptr_SC_PtrEgr_ortho_norm_UMI_matrix <-
  Ptr_SC_PtrEgr_ortho_UMI_matrix %>%
  apply(1, multiply_by, Ptr_SC_norm_factor) %>%
  Matrix() %>% t()
Egr_SC_PtrEgr_ortho_norm_UMI_matrix <-
  Egr_SC_PtrEgr_ortho_UMI_matrix %>%
  apply(1, multiply_by, Egr_SC_norm_factor) %>%
  Matrix() %>% t()

Ptr_SC_PtrEgr_mean_ortho_norm_UMI_matrix <-
  calculate_mean_counts(Ptr_SC_PtrEgr_ortho_norm_UMI_matrix,
                        Ptr_SC_cell_cluster_10)
Egr_SC_PtrEgr_mean_ortho_norm_UMI_matrix <-
  calculate_mean_counts(Egr_SC_PtrEgr_ortho_norm_UMI_matrix,
                        Egr_SC_cell_cluster_10)

output_as_10X_format(
  Ptr_SC_PtrEgr_ortho_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrEgr_ortho_UMI_matrix.rds"
)
output_as_10X_format(
  Egr_SC_PtrEgr_ortho_UMI_matrix,
  "output_10X_matrix/RDS_Egr_SC_PtrEgr_ortho_UMI_matrix.rds"
)
output_as_10X_format(
  Ptr_SC_PtrEgr_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrEgr_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Egr_SC_PtrEgr_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Egr_SC_PtrEgr_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Ptr_SC_PtrEgr_mean_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrEgr_mean_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Egr_SC_PtrEgr_mean_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Egr_SC_PtrEgr_mean_ortho_norm_UMI_matrix.rds"
)

## Construct the ortholog UMI matrix of Ptr and Tar
overlap_clusters <- intersect(unique(Ptr_ortho_table$clusterID),
                              unique(Tar_ortho_table$clusterID))
Ptr_SC_PtrTar_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Ptr_SC_UMI_matrix,
    Ptr_ortho_table,
    overlap_clusters
  )
Tar_SC_PtrTar_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Tar_SC_UMI_matrix,
    Tar_ortho_table,
    overlap_clusters
  )

Ptr_SC_PtrTar_ortho_norm_UMI_matrix <-
  Ptr_SC_PtrTar_ortho_UMI_matrix %>%
  apply(1, multiply_by, Ptr_SC_norm_factor) %>%
  Matrix() %>% t()
Tar_SC_PtrTar_ortho_norm_UMI_matrix <-
  Tar_SC_PtrTar_ortho_UMI_matrix %>%
  apply(1, multiply_by, Tar_SC_norm_factor) %>%
  Matrix() %>% t()

Ptr_SC_PtrTar_mean_ortho_norm_UMI_matrix <-
  calculate_mean_counts(Ptr_SC_PtrTar_ortho_norm_UMI_matrix,
                        Ptr_SC_cell_cluster_10)
Tar_SC_PtrTar_mean_ortho_norm_UMI_matrix <-
  calculate_mean_counts(Tar_SC_PtrTar_ortho_norm_UMI_matrix,
                        Tar_SC_cell_cluster_18)

output_as_10X_format(
  Ptr_SC_PtrTar_ortho_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrTar_ortho_UMI_matrix.rds"
)
output_as_10X_format(
  Tar_SC_PtrTar_ortho_UMI_matrix,
  "output_10X_matrix/RDS_Tar_SC_PtrTar_ortho_UMI_matrix.rds"
)
output_as_10X_format(
  Ptr_SC_PtrTar_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrTar_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Tar_SC_PtrTar_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Tar_SC_PtrTar_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Ptr_SC_PtrTar_mean_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrTar_mean_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Tar_SC_PtrTar_mean_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Tar_SC_PtrTar_mean_ortho_norm_UMI_matrix.rds"
)

## Construct the ortholog UMI matrix of Ptr and Lch
overlap_clusters <- intersect(unique(Ptr_ortho_table$clusterID),
                              unique(Lch_ortho_table$clusterID))
Ptr_SC_PtrLch_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Ptr_SC_UMI_matrix,
    Ptr_ortho_table,
    overlap_clusters
  )
Lch_SC_PtrLch_ortho_UMI_matrix <-
  calculate_ortholog_UMI_matrix(
    Lch_SC_UMI_matrix,
    Lch_ortho_table,
    overlap_clusters
  )

Ptr_SC_PtrLch_ortho_norm_UMI_matrix <-
  Ptr_SC_PtrLch_ortho_UMI_matrix %>%
  apply(1, multiply_by, Ptr_SC_norm_factor) %>%
  Matrix() %>% t()
Lch_SC_PtrLch_ortho_norm_UMI_matrix <-
  Lch_SC_PtrLch_ortho_UMI_matrix %>%
  apply(1, multiply_by, Lch_SC_norm_factor) %>%
  Matrix() %>% t()

Ptr_SC_PtrLch_mean_ortho_norm_UMI_matrix <-
  calculate_mean_counts(Ptr_SC_PtrLch_ortho_norm_UMI_matrix,
                        Ptr_SC_cell_cluster_10)
Lch_SC_PtrLch_mean_ortho_norm_UMI_matrix <-
  calculate_mean_counts(Lch_SC_PtrLch_ortho_norm_UMI_matrix,
                        Lch_SC_cell_cluster_10)

output_as_10X_format(
  Ptr_SC_PtrLch_ortho_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrLch_ortho_UMI_matrix.rds"
)
output_as_10X_format(
  Lch_SC_PtrLch_ortho_UMI_matrix,
  "output_10X_matrix/RDS_Lch_SC_PtrLch_ortho_UMI_matrix.rds"
)
output_as_10X_format(
  Ptr_SC_PtrLch_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrLch_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Lch_SC_PtrLch_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Lch_SC_PtrLch_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Ptr_SC_PtrLch_mean_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Ptr_SC_PtrLch_mean_ortho_norm_UMI_matrix.rds"
)
output_as_10X_format(
  Lch_SC_PtrLch_mean_ortho_norm_UMI_matrix,
  "output_10X_matrix/RDS_Lch_SC_PtrLch_mean_ortho_norm_UMI_matrix.rds"
)