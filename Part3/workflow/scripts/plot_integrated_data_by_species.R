library(Seurat)
library(scales)
library(magrittr)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
# Output figure
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

# Plot seurat UMAP colored on each sample
plot_on_sample <- function(
    MS_plotting_df,
    output_without_margin,
    output_without_legend,
    col_species,
    ...
) {
    x <- MS_plotting_df$UMAP.1
    y <- MS_plotting_df$UMAP.2
    col <- col_species[MS_plotting_df$Species]

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
            pch = 20, col = col_species, legend = names(col_species)
        )
    }
}


# Set parameters ===============================================================
# Get input parameters from command line
input_MS_plotting_csv <- snakemake@input$MS_plotting_csv
output_figure_folder <- snakemake@output$figure_folder


# Implementation ===============================================================
# Create output directory
if (!dir.exists(output_figure_folder)) {
    dir.create(output_figure_folder, recursive = TRUE)
}

# Input MS plotting information
MS_plotting_df <- read.csv(input_MS_plotting_csv)
n_species <- length(unique(MS_plotting_df$Species))

# Wrap function of plotting seurat UMAP colored on each sample
plotting_function <- function(
    output_without_margin,
    output_without_legend,
    col_species,
    ...
) {
    plot_on_sample(
        MS_plotting_df,
        output_without_margin,
        output_without_legend,
        col_species,
        ...
    )
}
## Multiples samples: colorful
if (n_species > 2) {
    col_matrix <- matrix(
        hue_pal()(n_species),
        ncol = n_species, nrow = n_species + 1,
        byrow = TRUE
    )
    for (i in seq(n_species)) {
        col_matrix[i, -i] <- paste0(col_matrix[i, -i], "00")
    }
}
## Two samples: (1st: black; 2nd: gold)
if (n_species == 2) {
    col_matrix <- matrix(
        c("#000000", "#C59739"),
        ncol = 2, nrow = 3,
        byrow = TRUE
    )
    for (i in seq(n_species)) {
        col_matrix[i, -i] <- paste0(col_matrix[i, i], "00")
    }
}
## One samples: black
if (n_species == 1) {
    col_matrix <- matrix("#000000", ncol = 1, nrow = 1)
}
## Output figure
for (i in seq(nrow(col_matrix))) {
    species <- unique(MS_plotting_df$Species)
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path =
            paste0(
                output_figure_folder,
                "/BySpecies_", c(species, "All")[i], ".png"
            ),
        output_without_margin = FALSE,
        output_without_legend = FALSE,
        col_species = col_matrix[i, ] %>%
            set_names(species)
    )
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path =
            paste0(
                output_figure_folder,
                "/BySpecies_", c(species, "All")[i], "_NoLegend.png"
            ),
        output_without_margin = FALSE,
        output_without_legend = TRUE,
        col_species = col_matrix[i, ] %>%
            set_names(species)
    )
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path =
            paste0(
                output_figure_folder,
                "/BySpecies_", c(species, "All")[i], "_Clear.png"
            ),
        output_without_margin = TRUE,
        output_without_legend = TRUE,
        col_species = col_matrix[i, ] %>%
            set_names(species)
    )
}
