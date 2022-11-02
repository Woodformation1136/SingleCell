library(magrittr)
library(getopt)
library(dplyr)


# Define functions =============================================================


# Set parameters ===============================================================
#' Get input parameters from command line
input_orthoMCL_txt <- snakemake@input$orthoMCL_txt
input_protein2gene_csv <- snakemake@input$protein2gene_csv
output_ortholog_long_csv <- snakemake@output$ortholog_long_csv
output_ortholog_long_convertedID_csv <-
    snakemake@output$ortholog_long_convertedID_csv


# Implementation ===============================================================
#' Create output directory
if (!dir.exists(dirname(output_ortholog_long_csv))) {
    dir.create(dirname(output_ortholog_long_csv), recursive = TRUE)
}

#' Create connection to input file
orthoMCL_con <- file(input_orthoMCL_txt, "r")
ortholog_long_df <- structure(
    .Data = matrix(nrow = 0, ncol = 3),
    dimnames = list(
        c(), c("cluster_id", "species_name", "gene_name")
    )
)
write.table(
    ortholog_long_df,
    output_ortholog_long_csv,
    row.names = FALSE, sep = ","
)

run_loop <- TRUE
run_loop_num <- 1
while (run_loop) {
    in_content <- readLines(orthoMCL_con, n = 1)

    run_loop <- (length(in_content) > 0)
    if (!run_loop) break

    temp <- in_content %>%
        gsub(pattern = '"', replacement = "") %>%
        strsplit(split = "\t") %>%
        extract2(1)
    cluster_id <- temp[1] %>%
        strsplit(split = "[(]") %>%
        extract2(1) %>%
        extract(1)
    gene_num <- temp[1] %>%
        strsplit(split = "[:|,]") %>%
        extract2(1) %>%
        extract(2) %>%
        as.numeric()
    temp_gene <- temp[-1] %>%
        strsplit(split = "[(]") %>%
        sapply(extract, 1)

    stopifnot(all(regexpr("_", temp_gene) == 4)) # TRUE
    species_name <- substr(temp_gene, 1, 3)
    gene_name <- substring(temp_gene, 5)

    ortholog_long_df <- data.frame(cluster_id, species_name, gene_name)
    # gene_num_vector[run_loop_num] <- gene_num

    write.table(
        ortholog_long_df,
        output_ortholog_long_csv,
        row.names = FALSE, col.names = FALSE,
        append = TRUE, quote = FALSE, sep = ","
    )

    cat("\r", run_loop_num)
    run_loop_num %<>% add(1)
}

close(orthoMCL_con)


#' Load the ortholog table (long data format)
ortho_df <- read.csv(output_ortholog_long_csv)

#' Convert the protein ID into gene ID for Ptr, Egr and Tar
protein2gene_df <- read.csv(input_protein2gene_csv)
ortho_convertedID_df <- ortho_df
ortho_convertedID_df$gene_name[
    match(protein2gene_df$protein, ortho_convertedID_df$gene_name)
] <- protein2gene_df$gene


write.csv(
    ortho_convertedID_df,
    file = output_ortholog_long_convertedID_csv,
    quote = FALSE,
    row.names = FALSE
)