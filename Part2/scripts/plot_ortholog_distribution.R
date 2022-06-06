library(magrittr)
library(tidyr)
library(dplyr)
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

# Get label encoding vector
# get_encoding_vector <- function(OrthoGroupNo, ortholog_long, species_order) {
#     select_species_name <- ortholog_long %>%
#         dplyr::filter(cluster_id == OrthoGroupNo) %>%
#         extract(, "species_name") %>%
#         unique()
#     out <- sapply(species_order, is_in, select_species_name)
#     out <- out %>% as.numeric() %>% set_names(species_order)
#     return(out)
# }

# Get label encoding vector
get_encoding_vector <- function(OrthoGroupNo, ortholog_long, species_order) {
    select_species_name <- ortholog_long %>%
        dplyr::filter(cluster_id %in% strsplit(OrthoGroupNo, "_")[[1]]) %>%
        magrittr::extract(, "species_name")
    out <-
        sapply(
            species_order,
            function(sp) return(sum(select_species_name == sp))
        ) %>%
        set_names(species_order)
    return(out)
}


# Set parameters ===============================================================
# Get input parameters from snakemake
ortholog_long_csv <- snakemake@input$ortholog_long_csv
plot_ortholog_csv <- snakemake@input$plot_ortholog_csv
output_dir <- snakemake@output$output_dir
species_order <- snakemake@params$species_order


# Implement program ============================================================
# Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input long ortholog table
ortholog_long <- read.csv(ortholog_long_csv)

# Input ortholog table for plotting
plot_ortholog <- read.csv(plot_ortholog_csv)
plot_cluster_id <- plot_ortholog$OrthoGroupNo
plot_x_axis_label <- plot_ortholog$CommonName

# Construct 0-or-1-matrix for plotting
plot_matrix <- sapply(
    plot_cluster_id,
    get_encoding_vector,
    ortholog_long = ortholog_long,
    species_order = species_order
)
n_ortholog <- ncol(plot_matrix)
n_species <- nrow(plot_matrix)

# # Construct metrics matrix for ordering
# ordering_matrics <- apply(
#     plot_matrix, MARGIN = 2,
#     function(ortho_vector) {
#         ortho_dist <- ifelse(ortho_vector > 0, 1, 0)
#         n_one <- sum(ortho_dist == 1)
#         i_first_one <- min(which(ortho_dist == 1), 14)
#         len_max_continuous_zero <-
#             max(nchar(unlist(strsplit(paste0(ortho_dist, collapse = ""), "1"))))
#         mean_n_ortho <- mean(ortho_vector)
#         out <- c(n_one, -i_first_one, len_max_continuous_zero, mean_n_ortho)
#         return(out)
#     }
# )

# # Reorder plotting data
# order_id <- order(
#     ordering_matrics[1, ],
#     ordering_matrics[2, ],
#     ordering_matrics[3, ],
#     ordering_matrics[4, ]
# )
# order_plot_cluste_id <- plot_cluster_id[order_id]
# order_plot_x_axis_label <- plot_x_axis_label[order_id]
# order_plot_matrix <- plot_matrix[, order_id]
order_plot_cluste_id <- plot_cluster_id
order_plot_x_axis_label <- plot_x_axis_label
order_plot_matrix <- plot_matrix

# Define plotting function for ortholog distribution
plotting_function <- function(output_with_legend, ...) {
    par(mar = c(9, 4, 1, 1) + 0.1)
    plot(NULL,
        xlim = c(0.5, n_ortholog + 0.5),
        ylim = c(0.5, n_species + 0.5),
        axes = F,
        panel.first = {
            box(bty = "L")
            axis(1, seq(n_ortholog), rep("", n_ortholog), las = 2)
            text(
                seq(n_ortholog), par("usr")[3] - 0.5,
                srt = 60, adj = 1, xpd = TRUE,
                labels = order_plot_x_axis_label
            )
            axis(2, seq(n_species), rev(species_order), las = 2)
        },
        xlab = "", # "Ortholog cluster id"
        ylab = "Species"
    )
    for (i in seq(n_ortholog)) {
        for (j in seq(14)) {
            n_ortholog <- order_plot_matrix[n_species + 1 - j, i]
            points(
                i, j,
                cex = ifelse(
                    test = n_ortholog > 0,
                    yes = log10(n_ortholog + 3),
                    no = 0 # 0.5
                ),
                pch = ifelse(
                    test = n_ortholog > 0,
                    yes = 19,
                    no = 1 # 4
                )
            )
        }
    }
    if (output_with_legend) {
        leg <- c(0, 1, 2, 5, 10, 20, 50)
        legend(
            "topright",
            xpd = TRUE, adj = c(0.5, 0.5),
            x.intersp = 2, y.intersp = 2,
            pch = ifelse(leg > 0, 19, 1),
            pt.cex = ifelse(leg > 0, log10(leg + 3), 0),
            legend = leg
        )
    }
    par(ori_par)
}

# Output figure
output_png_figure(
    plotting_function,
    output_figure = TRUE,
    output_path = paste0(
        output_dir, "/OrthologDistribution_Legend.png"
    ),
    output_without_margin = FALSE,
    output_with_legend = TRUE
)
output_png_figure(
    plotting_function,
    output_figure = TRUE,
    output_path = paste0(
        output_dir, "/OrthologDistribution.png"
    ),
    output_without_margin = FALSE,
    output_with_legend = FALSE
)
