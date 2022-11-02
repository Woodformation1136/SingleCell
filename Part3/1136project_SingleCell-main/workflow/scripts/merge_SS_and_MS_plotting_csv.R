suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tools)
})

# Define functions =============================================================


# Set parameters ===============================================================
#' Get input parameters from command line
input_MS_plotting_csv <- snakemake@input$MS_plotting_csv
input_SS_plotting_csvs <- snakemake@input$SS_plotting_csv_list
output_MS_plotting_csv <- snakemake@output$MS_plotting_csv


# Implementation ===============================================================
#' Input MS_plotting_df
MS_plotting_df <- read.csv(input_MS_plotting_csv)

#' Input SS_plotting_df list
SS_plotting_df_list <- sapply(
    simplify = FALSE,
    USE.NAMES = FALSE,
    input_SS_plotting_csvs,
    function(csv_path) {
        sample_name <- csv_path %>%
            basename() %>%
            file_path_sans_ext() %>%
            sub(pattern = "plotting_", replacement = "")
        out <- read.csv(csv_path, row.names = "Barcode")
        rownames(out) <- paste0(sample_name, "_", rownames(out))
        return(out)
    }
)
SS_plotting_df <- do.call(rbind, SS_plotting_df_list)
colnames(SS_plotting_df) <- paste0("SS_", colnames(SS_plotting_df))

#' Get sample_barcode in MS_plotting_df
sample_barcode <- paste0(MS_plotting_df$Sample, "_", MS_plotting_df$Barcode)

if (!all(sample_barcode %in% rownames(SS_plotting_df))) {
    stop("Some barcodes in MS_plotting_csv are not in SS_plotting_csv_list.")
}

#' Merge MS_plotting_df and SS_plotting_df
merge_plotting_df <- cbind(MS_plotting_df, SS_plotting_df[sample_barcode, ])

if (any(is.na(merge_plotting_df))) {
    stop("There is NA in merged plotting data frame.")
}

write.csv(
    merge_plotting_df,
    file = output_MS_plotting_csv,
    row.names = FALSE, quote = FALSE
)