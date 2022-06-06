library(magrittr)
library(tools)
library(rsvg)


# Define functions =============================================================
can_convert_to_numeric <- function(x) {
  all(grepl('^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$', x, perl = TRUE))
}
Add_single_color_tag <- function(text_line, sp_col, sp_name, marked_list) {
  text_list <- strsplit(text_line, "[<>]")[[1]]

  text_content <- gsub(pattern = " ", replacement = "", text_list[3])
  if (can_convert_to_numeric(text_content)) {
    new_content <- text_line
  } else {
    content_species <- substr(text_content, 1, 3)
    content_gene <- substring(text_content, 4)

    text_tail <- paste0("</text>")

    text_head <- paste0("<", text_list[2], ">")
    text_x <- strsplit(text_head, '"')[[1]][2] %>% as.numeric()
    text_x_left <- strsplit(text_head, '"')[[1]][1]
    text_x_right <-
      strsplit(text_head, '"')[[1]][-c(1, 2)] %>%
      paste(collapse = '"')

    if (content_gene %in% marked_list) {
      text_x_right <- sub("#000000", sp_col[content_species], text_x_right)
    }

    new_content <-
      paste0(
        sub("#000000", sp_col[content_species], text_head),
        paste0(" ", sp_name[content_species]),
        text_tail,
        paste(c(text_x_left, text_x + 35, text_x_right), collapse = '"'),
        content_gene,
        text_tail
      )
  }

  return(new_content)
}


# Set parameters ===============================================================
## Get input parameters from snakemake
expressed_protein_id_txt <- snakemake@input$expressed_protein_id_txt
input_tree_dir <- snakemake@input$input_tree_dir
output_dir <- snakemake@output$output_dir
# expressed_protein_id_txt <- "results/Expressed_proteinID/Multi_Species4_Sample5.txt"
# input_tree_dir <- "manual_results/Tree_from_MEGA"
# output_dir <- "results/Phylogenetic_tree_decorated/Multi_Species4_Sample5"

## Create species to color vector
sp_col <- c(
  PhP = "#000000", # "#D9D9D9",
  MaP = "#000000", # "#D9D9D9",
  SeM = "#000000", # "#595959",
  PiT = "#000000", # "#548235",
  GnM = "#000000", # "#548235",
  AmT = "#000000", # "#D65F47",
  OrS = "#000000", # "#D65F47",
  PoT = "#D65F47",
  EuG = "#D65F47",
  ArT = "#000000", # "#D65F47",
  CoC = "#000000", # "#D65F47",
  SoL = "#000000", # "#D65F47",
  LiC = "#FFC000",
  TrA = "#00B0F0"
)

## Conversion table for species name
sp_name <- c(
  PhP = "Ppa",
  MaP = "Mpo",
  SeM = "Smo",
  PiT = "Pta",
  GnM = "Gmo",
  AmT = "Atr",
  OrS = "Osa",
  PoT = "Ptr",
  EuG = "Egr",
  ArT = "Ath",
  CoC = "Cca",
  SoL = "Sly",
  LiC = "Lch",
  TrA = "Tar"
)


# Implementation ===============================================================
## Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

## Input expressed protein ID
expressed_protein_id <- readLines(expressed_protein_id_txt)

## Revise and save as new svg
file_list <-
  list.files(input_tree_dir) %>%
  grep(pattern = ".svg", value = TRUE) %>%
  grep(pattern = "_color", invert = TRUE, value = TRUE)

for (fin in file_list) {
  svg_lines <- readLines(paste0(input_tree_dir, "/", fin))
  svg_lines_revised <- sapply(
    svg_lines,
    function(x) {
      if (grepl("text", x)) {
        out <-
          Add_single_color_tag(
            text_line = x,
            sp_col = sp_col,
            sp_name = sp_name,
            marked_list = expressed_protein_id
          )
      } else {
        out <- x
      }
      return(out)
    }
  )
  output_filename <-
    paste0(output_dir, "/", file_path_sans_ext(fin), "_color")
  writeLines(
    svg_lines_revised,
    paste0(output_filename, ".svg")
  )
  rsvg_png(
    paste0(output_filename, ".svg"),
    paste0(output_filename, ".png")
  )
}
