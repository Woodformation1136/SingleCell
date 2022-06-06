library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(Matrix)
library(scales)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
# Output figure
output_png_figure <- function(
    plotting_function,
    # x, y, col, main = "",
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

# Plot seurat UMAP colored on each sample
plot_on_sample <- function(
    seurat_umap,
    output_without_margin,
    output_without_legend,
    col_sample,
    ...
) {
    x <- seurat_umap$UMAP_1
    y <- seurat_umap$UMAP_2
    col <- col_sample[seurat_umap$Sample]

    randam_order <- sample(length(x))
    plot(x[randam_order], y[randam_order],
        col = col[randam_order], pch = 20, cex = 0.3,
        xlab = "UMAP_1", ylab = "UMAP_2",
        main = ifelse(output_without_margin, "", "main"),
        axes = !output_without_margin, las = 1
    )
    if (!output_without_legend) {
        legend(
            "bottomleft",
            pch = 20, col = col_sample, legend = names(col_sample)
        )
    }
}


# Set parameters ===============================================================
# Get input parameters from command line
input_integrated_rds <- snakemake@input$input_integrated_rds
output_dir <- snakemake@output$output_dir


# Implementation ===============================================================
# Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input integrated data
combined_object <- readRDS(input_integrated_rds)
seurat_umap <- as.data.frame(
    combined_object@reductions$umap@cell.embeddings
)
seurat_umap$SampleBarcode <- rownames(seurat_umap)
seurat_umap$Sample <-
    seurat_umap$SampleBarcode %>%
    strsplit(split = "_") %>%
    sapply(extract, simplify = FALSE, 1:2) %>%
    sapply(paste, collapse = "_")
n_samples <- length(unique(seurat_umap$Sample))

# Wrap function of plotting seurat UMAP colored on each sample
plotting_function <- function(
    output_without_margin,
    output_without_legend,
    col_sample,
    ...
) {
    plot_on_sample(
        seurat_umap,
        output_without_margin,
        output_without_legend,
        col_sample,
        ...
    )
}
## Multiples samples: colorful
if (n_samples > 2) {
    col_matrix <- matrix(
        hue_pal()(n_samples), # c("#C59739", "#9B9B9B", "#000000"),
        ncol = n_samples, nrow = n_samples + 1,
        byrow = TRUE
    )
    for (i in seq(n_samples)) {
        col_matrix[i, -i] <- paste0(col_matrix[i, -i], "00")
    }
}
## Two samples: (1st: black; 2nd: gold)
if (n_samples == 2) {
    col_matrix <- matrix(
        c("#000000", "#C59739"),
        ncol = n_samples, nrow = n_samples + 1,
        byrow = TRUE
    )
    for (i in seq(n_samples)) {
        col_matrix[i, -i] <- paste0(col_matrix[i, i], "00")
    }
}
## Output figure
for (i in seq(nrow(col_matrix))) {
    sample_name <- unique(seurat_umap$Sample)
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path =
            paste0(output_dir, "/BySample_", c(sample_name, "All")[i], ".png"),
        output_without_margin = FALSE,
        output_without_legend = FALSE,
        col_sample = col_matrix[i, ] %>%
            set_names(sample_name)
    )
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path =
            paste0(
                output_dir,
                "/BySample_", c(sample_name, "All")[i], "_NoLegend.png"
            ),
        output_without_margin = FALSE,
        output_without_legend = TRUE,
        col_sample = col_matrix[i, ] %>%
            set_names(sample_name)
    )
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path =
            paste0(
                output_dir, "/BySample_", c(sample_name, "All")[i], "_Clear.png"
            ),
        output_without_margin = TRUE,
        output_without_legend = TRUE,
        col_sample = col_matrix[i, ] %>%
            set_names(sample_name)
    )
}
