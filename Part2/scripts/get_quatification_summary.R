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

get_statistics <- function(UMI_matrix) {
    out <- c(
        "n_cells" = sum(colSums(UMI_matrix) != 0),
        "minimum_UMI" = min(colSums(UMI_matrix)),
        "n_detected_genes" = sum(rowSums(UMI_matrix) != 0),
        "p_detected_genes(%)" =
            round(sum(rowSums(UMI_matrix) != 0) / nrow(UMI_matrix) * 100, 2),
        "median_genes_per_cell" = median(colSums(UMI_matrix != 0))
    )
    return(out)
}
## Mean reads per cell are required to be calculated manually


# Set parameters ===============================================================
## Get input parameters from snakemake
input_dirs <- snakemake@input$input_dirs
output_csv <- snakemake@output$output_csv
species_list <- snakemake@params$species_list

# input_dirs <- c(
#     "results/SCseq_TenX_Ptr",
#     "results/SCseq_TenX_Ptr2",
#     "results/SCseq_MARSseq_Ptr3",
#     "results/SCseq_MARSseq_Ptr4",
#     "results/SCseq_TenX_PalChen2021",
#     "results/SCseq_TenX_Lch",
#     "results/SCseq_MARSseq_Egr",
#     "results/SCseq_MARSseq_Egr2",
#     "results/SCseq_MARSseq_Tar",
#     "results/SCseq_MARSseq_Tar2",
#     "results/SCseq_TenX_Tar3"
# )
# species_list <- list(
#     "TenX_Ptr" = "Ptr_TenX",
#     "TenX_Ptr2" = "Ptr_TenX",
#     "MARSseq_Ptr3" = "Ptr_MARSseq",
#     "MARSseq_Ptr4" = "Ptr_MARSseq",
#     "TenX_PalChen2021" = "Pal",
#     "TenX_Lch" = "Lch",
#     "MARSseq_Egr" = "Egr",
#     "MARSseq_Egr2" = "Egr",
#     "MARSseq_Tar" = "Tar_MARSseq",
#     "MARSseq_Tar2" = "Tar_MARSseq",
#     "TenX_Tar3" = "Tar_TenX"
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
species_df <-
    data.frame(
        species = as.character(species_list),
        sample = names(species_list)
    ) %>%
    filter(sample %in% input_samples)
unique_species_vector <- unique(species_df$species)

# Merge matrices based on species
merge_UMI_matrix <- list()
for (sp in unique_species_vector) {
    sample_vector <- with(species_df, sample[species == sp])
    if (length(sample_vector) == 1) {
        merge_UMI_matrix[[sp]] <-
            input_UMI_matrix[[sample_vector]]
    } else {
        checkpoint <- sapply(input_UMI_matrix[sample_vector], rownames)
        stopifnot(is.matrix(checkpoint))
        ref_gene_id <- checkpoint[, 1]
        merge_UMI_matrix[[sp]] <-
            do.call(
                cbind,
                sapply(
                    input_UMI_matrix[sample_vector],
                    simplify = FALSE,
                    function(X) X[ref_gene_id, ]
                )
            )
    }
}

# Compute descriptive statistics of each species
descriptive_statistics <- sapply(merge_UMI_matrix, get_statistics)

# Output descriptive statistics as csv
write.csv(
    as.data.frame(descriptive_statistics),
    file = output_csv
)

# # Calculate mean reads per cell
# target_dir <- "rawdata/20220520_MAR-seq_organized/Reads_per_well/"
# child_dirs <- list.files(target_dir, full.names = TRUE)
# MARSseq_reads_file_list <-
#     sapply(child_dirs, function(x) paste0(x, "/", list.files(x))) %>%
#     unlist() %>%
#     setNames(NULL)

# MARSseq_well_reads <-
#     do.call(
#         rbind,
#         sapply(
#             MARSseq_reads_file_list,
#             simplify = FALSE,
#             function(fin) read.table(fin, header = TRUE)
#         )
#     ) %>%
#     set_rownames(NULL)
# dim(MARSseq_well_reads) # 30720 2

# ## Lch TenX (from result of cellranger)
# 357098221 / 2977 # 119952.4

# ## Ptr TenX (from result of cellranger)
# # 328948141 / 4705
# # 632895700 / 14302
# (328948141 + 632895700) / (4705 + 14302) # 50604.72

# ## Ptr MARS-seq (from result of cellranger)
# sum(with(MARSseq_well_reads, total[grepl("Pt", well_id)])) /
#     ncol(merge_UMI_matrix$Ptr_MARSseq) # 11837.35
# # 72906210 / 6159

# ## Egr MARS-seq
# sum(with(MARSseq_well_reads, total[grepl("Eg", well_id)])) /
#     ncol(merge_UMI_matrix$Egr) # 12910.85
# # 116313853 / 9009

# ## Tar MARS-seq
# sum(with(MARSseq_well_reads, total[grepl("Ta", well_id)])) /
#     ncol(merge_UMI_matrix$Tar_MARSseq) # 15406.75
# # 76602356 / 4972
