library(magrittr)
library(dplyr)
library(seqinr)


# Define functions =============================================================
get_species_peptide_id <- function(ortholog_id_string) {
    ortholog_id <- strsplit(ortholog_id_string, "_")[[1]]
    selected_df <- filter(ortholog_long, cluster_id %in% ortholog_id)
    out <- with(selected_df, paste(species_name, gene_name, sep = "_"))
    return(out)
}


# Set parameters ===============================================================
# Get input parameters from snakemake
combined_faa <- snakemake@input$combined_faa
ortholog_long_csv <- snakemake@input$ortholog_long_csv
plot_ortholog_csv <- snakemake@params$plot_ortholog_csv
output_dir <- snakemake@output$output_dir



# Implement program ============================================================
# Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input data
ortholog_long <- read.csv(ortholog_long_csv)
combined_fasta <- read.fasta(
    combined_faa, seqtype = "AA", as.string = TRUE
)
plot_ortholog <- read.csv(plot_ortholog_csv)
combined_name <- names(combined_fasta)

# Get and write fasta with ortholog id
for (i in seq(nrow(plot_ortholog))) {
    selected_peptide_id <-
        get_species_peptide_id(plot_ortholog$OrthoGroupNo[i])
    stopifnot(all(selected_peptide_id %in% combined_name))

    selected_index <- match(selected_peptide_id, combined_name)
    write.fasta(
        sequences = combined_fasta[selected_index],
        names = names(combined_fasta[selected_index]),
        as.string = TRUE,
        file.out = paste0(
            output_dir, "/",
            gsub("*", "", plot_ortholog$CommonName[i], fixed = TRUE),
            ".faa"
        )
    )
}

