library(magrittr)
library(tidyr)
library(dplyr)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================


# Set parameters ===============================================================
# Get input parameters from snakemake
all_ortholog_long_csv <- snakemake@input[[1]]
primary_ortholog_long_csv <- snakemake@output[[1]]
gff_folder <- file.path(
    snakemake@params$rawdata_path,
    "20220000_Single_cell_rawdata/gff_to_get_one_peptide_per_gene"
)


# Implement program ============================================================
# Input long ortholog table
all_ortholog_long <- read.csv(all_ortholog_long_csv)

# Remove redundant peptides (to one peptide per gene)
## 1. Get transcript/CDS id with longest CDS for each gene
## Arabidopsis thaliana
gff_ArT <- read.table(
    paste0(gff_folder, "/Araport11_GFF3_genes_transposons.201606.gff"),
    sep = "\t", quote = ""
)
transcript2gene <- gff_ArT %>%
    filter(V3 == "mRNA") %>%
    select(V9) %>%
    unique() %>%
    separate(
        col = V9, into = c("transcript", "gene"),
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        transcript = sub("ID=", "", transcript),
        gene = sub("Parent=", "", gene)
    )
CDS2transcript <- gff_ArT %>%
    filter(V3 == "CDS") %>%
    separate(
        col = V9, into = c("CDS", "transcript"),
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        transcript = sub("Parent=", "", transcript)
    ) %>%
    group_by(transcript) %>%
    group_modify(function(.x, .y) {
        return(data.frame(len_CDS = sum(.x$V5 - .x$V4)))
    })
CDS2gene <-
    full_join(CDS2transcript, transcript2gene, by = "transcript") %>%
    na.omit()
all_ortholog_long %>%
    filter(species_name == "ArT") %>%
    magrittr::extract(, "gene_name") %>%
    setequal(CDS2gene$transcript) %>%
    stopifnot()
keep_id_ArT <- CDS2gene %>%
    group_by(gene) %>%
    group_map(function(.x, .y) {
        return(with(.x, transcript[which.max(len_CDS)]))
    }) %>%
    unlist()
## Oryza sativa
gff_OrS <- read.table(
    paste0(gff_folder, "/Oryza_sativa.IRGSP-1.0.52.gff3"),
    sep = "\t", quote = ""
)
transcript2gene <- gff_OrS %>%
    filter(V3 == "mRNA") %>%
    select(V9) %>%
    unique() %>%
    separate(
        col = V9, into = c("transcript", "gene"),
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        transcript = sub("ID=transcript:", "", transcript),
        gene = sub("Parent=gene:", "", gene)
    )
CDS2transcript <- gff_OrS %>%
    filter(V3 == "CDS") %>%
    separate(
        col = V9, into = c("CDS", "transcript"),
        sep = ";", extra = "drop"
    ) %>%
    mutate(
        CDS = sub("ID=CDS:", "", CDS),
        transcript = sub("Parent=transcript:", "", transcript)
    ) %>%
    group_by(transcript) %>%
    group_modify(function(.x, .y) {
        return(data.frame(CDS = .x$CDS[1], len_CDS = sum(.x$V5 - .x$V4)))
    })
CDS2gene <-
    full_join(CDS2transcript, transcript2gene, by = "transcript") %>%
    na.omit()
all_ortholog_long %>%
    filter(species_name == "OrS") %>%
    magrittr::extract(, "gene_name") %>%
    setequal(CDS2gene$CDS) %>%
    stopifnot()
keep_id_OrS <- CDS2gene %>%
    group_by(gene) %>%
    group_map(function(.x, .y) {
        return(with(.x, CDS[which.max(len_CDS)]))
    }) %>%
    unlist()
## Marchantia polymorpha
gff_MaP <- read.table(
    paste0(gff_folder, "/MpTak1v5.1_r1.gtf"),
    sep = "\t", quote = ""
)
CDS2gene <- gff_MaP %>%
    filter(V3 == "CDS") %>%
    separate(
        col = V9, into = c("transcript", "gene"),
        sep = "\";", extra = "drop"
    ) %>%
    mutate(
        transcript = sub("transcript_id \"", "", transcript),
        gene = sub(" gene_id \"", "", gene)
    ) %>%
    group_by(transcript) %>%
    group_modify(function(.x, .y) {
        return(data.frame(gene = .x$gene[1], len_CDS = sum(.x$V5 - .x$V4)))
    })
all_ortholog_long %>%
    filter(species_name == "MaP") %>%
    magrittr::extract(, "gene_name") %>%
    setequal(CDS2gene$transcript) %>%
    stopifnot()
keep_id_MaP <- CDS2gene %>%
    group_by(gene) %>%
    group_map(function(.x, .y) {
        return(with(.x, transcript[which.max(len_CDS)]))
    }) %>%
    unlist()
## 2. Keep only selected ids in ortholog table
primary_ortholog_long <-
    all_ortholog_long %>%
    magrittr::extract(
        -which(with(
            all_ortholog_long,
            species_name %in% c("ArT", "OrS", "MaP") &
                !(gene_name %in% c(keep_id_ArT, keep_id_OrS, keep_id_MaP))
        )),
    )
# dim(primary_ortholog_long) # 478525 3 -> 445911 3


# Output primary ortholog table
write.csv(
    primary_ortholog_long,
    file = primary_ortholog_long_csv,
    quote = FALSE
)