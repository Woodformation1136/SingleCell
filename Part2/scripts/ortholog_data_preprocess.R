library(magrittr)
library(getopt)
library(dplyr)

spec <- matrix(
  c(
    "input", "i", 2, "character", "Input file",
    "output", "o", 2, "character", "Output file",
    "output_convertedID", "t", 2, "character", "Output file after ID convertion"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)
input <- opt$input
output <- opt$output
output_convertedID <- opt$output_convertedID

in_file <- file(input, "r")
out_file <- structure(
  .Data = matrix(nrow = 0, ncol = 3),
  dimnames = list(c(), c("cluster_id", "species_name", "gene_name"))
)
write.table(out_file, output,
  row.names = FALSE, sep = ","
)

# gene_num_vector <- c()

run_loop <- TRUE
run_loop_num <- 1
while (run_loop) {
  in_content <- readLines(in_file, n = 1)

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

  out_dataframe <- data.frame(cluster_id, species_name, gene_name)
  # gene_num_vector[run_loop_num] <- gene_num

  write.table(
    out_dataframe, output,
    row.names = FALSE, col.names = FALSE,
    append = TRUE, quote = FALSE, sep = ","
  )

  cat("\r", run_loop_num)
  run_loop_num %<>% add(1)
}

close(in_file)

# sum(gene_num_vector)


## Load the ortholog table (long data format)
ortho_table <- read.csv(output)
# dim(ortho_table) # 478525 3

## Convert the protein ID into gene ID
### Define the function for ID extracting
extract_prefix_before_point <- function(raw_id, which_point = 2) {
  cut_position <-
    raw_id %>%
    gregexpr(pattern = "[.]") %>%
    extract2(1) %>%
    extract(which_point)
  extracted_id <- substr(raw_id, 1, cut_position - 1)
  extracted_id
}
### Replace Ptr peptide ID with gene ID
Ptr_protein_fasta <-
  scan(paste0(
    "/home/woodformation/HDD1/GenomicsData/Ptrichocarpa_533_v4.1/",
    "Ptrichocarpa_533_v4.1.protein_primaryTranscriptOnly.fa"
  ),
  what = "", sep = "\n"
  )
fasta_info <- grep(">", Ptr_protein_fasta, value = TRUE)
Ptr_protein_ID <-
  strsplit(fasta_info, " ") %>%
  sapply(extract, 1) %>%
  sub(">", "", .)
Ptr_gene_ID <-
  strsplit(fasta_info, " ") %>%
  sapply(extract, 5) %>%
  sub("ID=", "", .) %>%
  sapply(extract_prefix_before_point, 2)
stopifnot(
  all(filter(ortho_table, species_name == "PoT")$gene_name %in% Ptr_protein_ID)
)
Ptr_row_indices <- which(ortho_table$species_name == "PoT")
ortho_table[Ptr_row_indices, "gene_name"] <-
  paste0(
    Ptr_gene_ID[match(
      ortho_table[Ptr_row_indices, "gene_name"],
      Ptr_protein_ID
    )],
    ".v4.1"
  )

### Replace Egr peptide ID with gene ID
Egr_row_indices <- which(ortho_table$species_name == "EuG")
ortho_table[Egr_row_indices, "gene_name"] <-
  sub("1.p", "v2.0", ortho_table[Egr_row_indices, "gene_name"])

### Replace Tar peptide ID with gene ID
Tar_gtf <-
  read.table(paste0(
    "/home/woodformation/HDD1/GenomicsData/Taralioides_20200702/",
    "Trochodendron_aralioides_chromosomes_pasa2.longest.filter.gtf"
  ), sep = "\t", header = FALSE)
ID_info <- unique(Tar_gtf$V9)
Tar_transcript_id <-
  strsplit(ID_info, ";") %>%
  sapply(extract, 1) %>%
  sub("transcript_id ", "", .)
Tar_gene_id <-
  strsplit(ID_info, ";") %>%
  sapply(extract, 2) %>%
  sub(" gene_id ", "", .)
  # sub("evm.TU.", "", .)
  # Dr. Ku remove prefix of gene id during quantification
stopifnot(
  all(
    filter(ortho_table, species_name == "TrA")$gene_name %in% Tar_transcript_id
  )
)
Tar_row_indices <- which(ortho_table$species_name == "TrA")
ortho_table[Tar_row_indices, "gene_name"] <-
  Tar_gene_id[match(
    ortho_table[Tar_row_indices, "gene_name"],
    Tar_transcript_id
  )]

write.csv(
  ortho_table,
  file = output_convertedID,
  quote = FALSE,
  row.names = FALSE
)