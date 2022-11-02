suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tools)
    library(RColorBrewer)
    ori_par <- par(no.readonly = TRUE)
})


# Define functions =============================================================
#' Output figure
output_png_figure <- function(
    plotting_function,
    output_figure = FALSE,
    output_path = "temp.png",
    output_without_margin = FALSE,
    ...
) {
    if (output_figure) {
        png(output_path,
            pointsize = 10, res = 300,
            width = 20, height = 15, units = "cm")
    }

    if (output_without_margin) {
        par(mai = c(0, 0, 0, 0))
    } else {
        par(mai = ori_par$mai)
    }

    plotting_function(
        output_figure = output_figure,
        output_path = output_path,
        output_without_margin = output_without_margin,
        ...
    )

    par(mai = ori_par$mai)

    if (output_figure) dev.off()
}

#' Plot seurat UMAP colored by correlation
plot_with_correlation <- function(
    MS_plotting_df = merge_plotting_df,
    cor_colname,
    output_without_margin,
    ...
) {
    x <- MS_plotting_df$UMAP.1
    y <- MS_plotting_df$UMAP.2
    cor_vector <- MS_plotting_df[, cor_colname]

    upper_bound <- quantile(cor_vector, 0.95)
    lower_bound <- upper_bound * 0.8
    color_index <-
        ((cor_vector - lower_bound) / (upper_bound - lower_bound)) %>%
        pmax(0) %>%
        pmin(1) %>%
        multiply_by(500) %>%
        round(digits = 0) %>%
        add(1)

    color_tick <- rev(brewer.pal(11, "Spectral"))
    color_pool <- colorRampPalette(color_tick[7:11])(501)

    randam_order <- sample(length(x))
    plot(x[randam_order], y[randam_order],
        col = color_pool[color_index[randam_order]],
        pch = 20, cex = 0.5,
        xlab = "UMAP_1", ylab = "UMAP_2",
        main = ifelse(
            output_without_margin, "",
            paste0(
                "Min:", round(min(cor_vector), 2), " ",
                "LB:", round(lower_bound, 2), " ",
                "UB:", round(upper_bound, 2), " ",
                "Max:", round(max(cor_vector), 2)
            )
        ),
        axes = !output_without_margin, las = 1
    )
}


# Set parameters ===============================================================
#' Get input parameters from command line
input_MS_plotting_csv <- snakemake@input$MS_plotting_csv
input_SS_MQD_plotting_csv <- snakemake@input$SS_MQD_plotting_csv
SC_sample <- snakemake@wildcards$SC_sample
output_figure_folder <- snakemake@output$figure_folder


# Implementation ===============================================================
#' Create output directory
if (!dir.exists(output_figure_folder)) {
    dir.create(output_figure_folder, recursive = TRUE)
}

#' Input plotting information
MS_plotting_df <- read.csv(input_MS_plotting_csv)
SS_MQD_plotting_df <- read.csv(input_SS_MQD_plotting_csv)
is_barcode_column <- grepl("Barcode", colnames(SS_MQD_plotting_df))
is_correlation_column <- grepl("Correlation_", colnames(SS_MQD_plotting_df))

#' Merge correlation into plotting data frame
selected_MS_plotting_df <- filter(MS_plotting_df, Sample == SC_sample)
if (nrow(selected_MS_plotting_df) == 0) {
    stop("No barcodes in SC_sample are in MS_plotting_csv.")
}
merge_plotting_df <- full_join(
    selected_MS_plotting_df,
    SS_MQD_plotting_df[, is_barcode_column | is_correlation_column],
    by = "Barcode"
)

# Output figures
correlation_columns <- colnames(SS_MQD_plotting_df)[is_correlation_column]
for (selected_column in correlation_columns) {
    output_png_figure(
        plotting_function = plot_with_correlation,
        MS_plotting_df = merge_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, ".png"),
        output_without_margin = FALSE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        MS_plotting_df = merge_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Clear.png"),
        output_without_margin = TRUE
    )
}
