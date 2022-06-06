library(tidyr)
library(dplyr)
library(magrittr)


# Define functions =============================================================
## Define the function for ID extracting
extract_prefix_before_point <- function(raw_id, which_point = 2) {
    cut_position <-
        raw_id %>%
        gregexpr(pattern = "[.]") %>%
        extract2(1) %>%
        extract(which_point)
    extracted_id <- substr(raw_id, 1, cut_position - 1)
    extracted_id
}


# Set parameters ===============================================================
## Get input parameters from snakemake
protein2gene_csv <- snakemake@output[[1]]


# Implementation ===============================================================
## Ptr protein ID to gene ID
fasta_protein_Ptr <-
    scan(
        paste0(
            "/home/woodformation/HDD1/GenomicsData/Ptrichocarpa_533_v4.1/",
            "Ptrichocarpa_533_v4.1.protein_primaryTranscriptOnly.fa"
        ),
        what = "", sep = "\n"
    )
fasta_info <- grep(">", fasta_protein_Ptr, value = TRUE)
protein_id_Ptr <-
    strsplit(fasta_info, " ") %>%
    sapply(extract, 1) %>%
    sub(">", "", .)
gene_id_Ptr <-
    strsplit(fasta_info, " ") %>%
    sapply(extract, 5) %>%
    sub("ID=", "", .) %>%
    sapply(extract_prefix_before_point, 2) %>%
    paste0(".v4.1")
protein2gene_Ptr <-
    data.frame(
        protein = protein_id_Ptr,
        gene = gene_id_Ptr
    )

## Egr protein ID to gene ID
gff_Egr <-
    read.table(
        paste0(
            "/home/woodformation/HDD1/GenomicsData/Egrandis_297_v2.0/",
            "Egrandis_297_v2.0.gene.gff3"
        ),
        sep = "\t", header = FALSE
    )
protein2gene_Egr <- gff_Egr %>%
    filter(V3 == "gene") %>%
    select(V9) %>%
    unique() %>%
    separate(
        col = V9, into = "transcript",
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        transcript = sub("ID=", "", transcript)
    ) %>%
    mutate(
        protein = sub("v2.0", "1.p", transcript),
        gene = transcript,
        .keep = "none"
    )

## Tar protein ID to gene ID
gtf_Tar <-
    read.table(
        paste0(
            "/home/woodformation/HDD1/GenomicsData/Taralioides_20200702/",
            "Trochodendron_aralioides_chromosomes_pasa2.longest.filter.gtf"
        ),
        sep = "\t", header = FALSE
    )
protein2gene_Tar <- gtf_Tar %>%
    filter(V3 == "transcript") %>%
    select(V9) %>%
    unique() %>%
    separate(
        col = V9, into = c("protein", "transcript"),
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        protein = sub("transcript_id ", "", protein),
        gene =
            transcript %>%
            sub(" gene_id ", "", .) %>%
            sub("evm.TU.", "", .),
            # Dr. Ku remove prefix of gene id during quantification
        .keep = "none"
    )

## Lch protein ID to gene ID
gff_Lch <-
    read.table(
        paste0(
            "/home/woodformation/HDD1/GenomicsData/",
            "Liriodendron_chinense_20191212/",
            "Liriodendron_chinense_gff"
        ),
        sep = "\t", header = FALSE
    )
protein2gene_Lch <- gff_Lch %>%
    filter(V3 == "mRNA") %>%
    select(V9) %>%
    unique() %>%
    separate(
        col = V9, into = "transcript",
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        protein = sub("ID=", "", transcript),
        gene = sub("ID=", "", transcript),
        .keep = "none"
    )

## Output protein ID to gene ID data frame
df_list <-
    list(
        protein2gene_Ptr,
        protein2gene_Egr,
        protein2gene_Tar,
        protein2gene_Lch
    )
write.csv(
    cbind(
        Species = rep(
            c("PoT", "EuG", "TrA", "LiC"),
            sapply(df_list, nrow)
        ),
        do.call(rbind, df_list)
    ),
    file = protein2gene_csv,
    row.names = FALSE
)