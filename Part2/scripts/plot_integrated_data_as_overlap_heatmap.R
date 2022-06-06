library(igraph)
library(magrittr)
library(Matrix)
library(MASS)
library(ggplot2)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
# Output figure
output_png_figure <- function(
    plotting_function,
    # x, y, col, main = "",
    output_figure = FALSE,
    output_path = "temp.png",
    output_clear = FALSE,
    output_ggplot = FALSE,
    ...
) {
    if (!output_ggplot) {
        if (output_figure) {
            png(output_path,
                pointsize = 10, res = 300,
                width = 20, height = 15, units = "cm"
            )
        }

        if (output_clear) {
            par(mai = c(0, 0, 0, 0))
        } else {
            par(mai = ori_par$mai)
        }

        plotting_function(
            output_figure = output_figure,
            output_path = output_path,
            output_clear = output_clear,
            ...
        )

        par(mai = ori_par$mai)

        if (output_figure) dev.off()
    } else {
        plotting_function(
            output_figure = output_figure,
            output_path = output_path,
            output_clear = output_clear,
            ...
        )

        if (output_figure) {
            ggsave(
                output_path,
                dpi = 300,
                width = 20, height = 15, units = "cm"
            )
        }
    }
}

# Contruct MST tree, cut inter-species edges, and extract center of each tree
get_mst_subtree_center <- function(
    projection,
    ref_sample = "TenX_Ptr",
    all_sample
) {
    message("Calculate the distance between each pair of cells")
    dist_matrix <- as.matrix(dist(projection)) %>% Matrix(sparse = TRUE)
    stopifnot(sum(dist_matrix == 0) == nrow(projection))

    message("Create the graph from adjacent matrix")
    graph_full <- graph_from_adjacency_matrix(
        dist_matrix,
        mode = "undirected",
        weighted = TRUE
    )

    message("Construct the MST")
    graph_mst <- mst(graph_full)

    message("Remove inter-species edges")
    edge_vname <- attr(E(graph_mst), "vnames")
    # is_cross_sample_edge <- edge_vname %>%
    #     gregexpr(pattern = paste0(ref_sample, "_")) %>%
    #     sapply(function(x) (length(x) == 1) & (x[1] != -1))
    vertex_matrix <- sapply(
        paste0(all_sample, "_"),
        function(pattern_sample) {
            out <- unlist(regexec(
                pattern = pattern_sample,
                text = edge_vname
            ))
            return(out != -1)
        }
    )
    stopifnot(all(rowSums(vertex_matrix) %in% 1:2))
    is_cross_sample_edge <- (rowSums(vertex_matrix) == 2)
    del_edge <- edge_vname[is_cross_sample_edge]
    graph_cut_mst <- delete_edges(graph_mst, del_edge)

    message("Extract the subgraph centers")
    candidate_vertices <- attr(V(graph_cut_mst), "name")
    subgraph_center <- rep(NA, length(candidate_vertices))
    i <- 1
    while (length(candidate_vertices) > 0) {
        picked_vertex <- candidate_vertices[1]
        picked_vertices <- attr(
            subcomponent(graph_cut_mst, picked_vertex),
            which = "name"
        )
        picked_graph <- induced_subgraph(graph_cut_mst, picked_vertices)
        picked_closeness <- closeness(picked_graph)
        if (length(picked_closeness) == 1) {
            subgraph_center[i] <- names(picked_closeness)
        } else {
            picked_center <- which(picked_closeness == max(picked_closeness))
            subgraph_center[i] <- sample(names(picked_center), 1)
        }
        candidate_vertices <- setdiff(candidate_vertices, picked_vertices)
        i <- i + 1
    }
    subgraph_center <- na.omit(subgraph_center)
    return(subgraph_center)
}

# Plot contour plot with center points
center_contour_plot <- function(
    subgraph_center = subtree_center,
    projection_umap = seurat_umap[, 1:2],
    ref_sample = "TenX_Ptr",
    lower_percentile = 5,
    higher_percentile = 40,
    main = ""
) {
    ref_projection_umap <-
        projection_umap[
            grepl(paste0(ref_sample, "_"), rownames(projection_umap)),
        ]
    ref_density_map <- kde2d(
        ref_projection_umap[, 1],
        ref_projection_umap[, 2],
        n = 500,
        lims = c(
            min(projection_umap[, 1]) - 0.5,
            max(projection_umap[, 1]) + 0.5,
            min(projection_umap[, 2]) - 0.5,
            max(projection_umap[, 2]) + 0.5
        )
    )
    get_territory_density <- function(coordinate) {
        i_x1 <- max(which(ref_density_map$x < coordinate[1]))
        i_y1 <- max(which(ref_density_map$y < coordinate[2]))
        df <- data.frame(
            X = ref_density_map$x[i_x1 + c(0, 0, 1, 1)],
            Y = ref_density_map$y[i_y1 + c(0, 1, 0, 1)],
            Z = apply(
                cbind(i_x1 + c(0, 0, 1, 1), i_y1 + c(0, 1, 0, 1)), 1,
                function(xy) ref_density_map$z[xy[1], xy[2]]
            )
        )
        out <- predict.lm(
            lm(Z ~ X + Y, data = df),
            data.frame(X = coordinate[1], Y = coordinate[2])
        )
        # out <- ref_density_map$z[
        #     max(which(ref_density_map$x < coordinate[1])),
        #     max(which(ref_density_map$y < coordinate[2]))
        # ]
        return(out)
    }
    ref_territory_density <- ref_projection_umap %>%
        apply(1, get_territory_density)
    norm_factor <- 1 / sum(ref_territory_density)

    ref_subgraph_center <- subgraph_center %>%
        extract(grepl(paste0(ref_sample, "_"), .))
    oth_subgraph_center <- subgraph_center %>%
        extract(!grepl(paste0(ref_sample, "_"), .))

    n_ref <- sum(grepl(paste0(ref_sample, "_"), rownames(projection_umap)))
    n_oth <- sum(!grepl(paste0(ref_sample, "_"), rownames(projection_umap)))
    n_total <- nrow(projection_umap)
    n_center_ref <- length(ref_subgraph_center)
    n_center_oth <- length(oth_subgraph_center)

    get_subgraph_center_density <- function(partsubgraph_center) {
        all_center_coordinate <- projection_umap[subgraph_center, ]
        part_center_coordinate <- projection_umap[partsubgraph_center, ]
        density_map <- kde2d(
            part_center_coordinate[, 1],
            part_center_coordinate[, 2],
            n = 500,
            h = apply(all_center_coordinate, 2, bandwidth.nrd) / 2,
            lims = c(
                min(projection_umap[, 1]) - 0.5,
                max(projection_umap[, 1]) + 0.5,
                min(projection_umap[, 2]) - 0.5,
                max(projection_umap[, 2]) + 0.5
            )
        )
        return(density_map)
    }
    ref_subgraph_center_density <-
        get_subgraph_center_density(ref_subgraph_center)
    oth_subgraph_center_density <-
        get_subgraph_center_density(oth_subgraph_center)
    stopifnot(ref_subgraph_center_density$x == oth_subgraph_center_density$x)
    stopifnot(ref_subgraph_center_density$y == oth_subgraph_center_density$y)

    merge_subgraph_center_density <- ref_subgraph_center_density
    merge_subgraph_center_density$z <-
        (n_ref / n_total) * n_center_ref * ref_subgraph_center_density$z +
        (n_oth / n_total) * n_center_oth * oth_subgraph_center_density$z
    merge_subgraph_center_density$z <-
        merge_subgraph_center_density$z %>% multiply_by(norm_factor)

    message(
        "Plot density max:",
        round(max(merge_subgraph_center_density$z), 5)
    )
    message(
        "Plot density Q75:",
        round(quantile(merge_subgraph_center_density$z, 0.75), 5)
    )
    message(
        "Plot density Q50:",
        round(quantile(merge_subgraph_center_density$z, 0.50), 5)
    )
    message(
        "Plot density Q25:",
        round(quantile(merge_subgraph_center_density$z, 0.25), 5)
    )
    message(
        "Plot density min:",
        round(min(merge_subgraph_center_density$z), 5)
    )

    territory_map <- kde2d(
        projection_umap[, 1],
        projection_umap[, 2],
        n = 500, h = 0.02,
        lims = c(
            min(projection_umap[, 1]) - 0.5,
            max(projection_umap[, 1]) + 0.5,
            min(projection_umap[, 2]) - 0.5,
            max(projection_umap[, 2]) + 0.5
        )
    )

    plot(NA,
        xlim = range(projection_umap[, "UMAP_1"]),
        ylim = range(projection_umap[, "UMAP_2"]),
        xlab = "", ylab = "", axes = FALSE, main = main
    )
    box()
    n_col_levels <- 500 # number of levels between 0 and 1
    lower_rank <- n_col_levels * lower_percentile / 100
    higher_rank <- n_col_levels * higher_percentile / 100
    if (max(merge_subgraph_center_density$z) > 1) {
        n_larger_than_1 <-
            (max(merge_subgraph_center_density$z) - 1) * n_col_levels
        filled_col <- c(
            colorRampPalette(
                c("#58b94e", "#fec84a")
            )(lower_rank),
            colorRampPalette(
                c("#fec84a", "#eb473e")
            )(higher_rank - lower_rank),
            colorRampPalette(
                c("#eb473e", "#713020")
            )(n_col_levels - higher_rank),
            rep("#713020", n_larger_than_1)
        )
        .filled.contour(
            merge_subgraph_center_density$x,
            merge_subgraph_center_density$y,
            merge_subgraph_center_density$z,
            levels = c(0, seq(n_col_levels + n_larger_than_1) / n_col_levels),
            col = filled_col
        )
        col_counts <-
            cut(
                merge_subgraph_center_density$z[territory_map$z > 0],
                breaks =
                    c(0, seq(n_col_levels + n_larger_than_1) / n_col_levels)
            ) %>%
            table() %>%
            as.vector()
        col_counts_df <- data.frame(filled_col = filled_col, col_counts)
    } else {
        filled_col <- c(
            colorRampPalette(
                c("#58b94e", "#fec84a")
            )(lower_rank),
            colorRampPalette(
                c("#fec84a", "#eb473e")
            )(higher_rank - lower_rank),
            colorRampPalette(
                c("#eb473e", "#713020")
            )(n_col_levels - higher_rank)
        )
        .filled.contour(
            merge_subgraph_center_density$x,
            merge_subgraph_center_density$y,
            merge_subgraph_center_density$z,
            levels = seq(0, 1, length.out = n_col_levels + 1),
            col = filled_col
        )
        col_counts <-
            cut(
                merge_subgraph_center_density$z[territory_map$z > 0],
                breaks = seq(0, 1, length.out = n_col_levels + 1)
            ) %>%
            table() %>%
            as.vector()
        col_counts_df <- data.frame(filled_col, col_counts)
    }

    .filled.contour(
        territory_map$x,
        territory_map$y,
        ifelse(territory_map$z > 0, 1, 0),
        levels = c(0, 0.5), col = c("white", NA)
    )
    contour(
        territory_map$x,
        territory_map$y,
        ifelse(territory_map$z > 0, 1, 0),
        levels = 0.5, lwd = 0.8, drawlabels = FALSE, add = TRUE
    )

    return(col_counts_df)
    # COL = rep("gray90",nrow(projection_umap))
    # COL[match(subgraph_center,rownames(projection_umap))] = 2
    # colInd = order(COL,decreasing=T)
    # points(projection_umap[colInd,"UMAP_1"],
    #        projection_umap[colInd,"UMAP_2"],
    #        col=COL[colInd],pch=20,cex=1)
}

# Plot distribution of col counts
plot_col_counts <- function(col_counts_df, type = "hist", df_path = NULL) {
    n_col_levels <- nrow(col_counts_df)
    agg_filled_col <- sapply(
        seq(n_col_levels / 5),
        function(i) {
            out <- col_counts_df$filled_col[min((i - 1) * 5 + 3, n_col_levels)]
            return(out)
        }
    )
    agg_col_counts <- sapply(
        seq(n_col_levels / 5),
        function(i) {
            out <- sum(
                col_counts_df$col_counts[
                    seq(1 + (i - 1) * 5, min(5 + (i - 1) * 5, n_col_levels))
                ]
            )
            return(out)
        }
    )
    col_counts_df <- data.frame(
        filled_col = agg_filled_col,
        col_counts = agg_col_counts
    )
    col_counts_df$col_props <-
        col_counts_df$col_counts / sum(col_counts_df$col_counts)

    if (!is.null(df_path)) write.csv(col_counts_df, df_path, row.names = FALSE)

    if (type == "hist") {
        plot(
            col_counts_df$col_props,
            col = col_counts_df$filled_col,
            type = "h", axes = FALSE,
            xlab = "Color",
            ylab = "Proportion (%)",
            cex.lab = 1.4,
            lwd = 3
        )
        box(bty = "L", lwd = 2)
        axis(
            side = 2,
            at = pretty(seq(0, max(col_counts_df$col_props) * 1000)) / 1000,
            labels = pretty(seq(0, max(col_counts_df$col_props) * 1000)) / 10,
            las = 2
        )
    } else if (type == "pie") {
        factor_filled_col <-
            with(col_counts_df, paste0(filled_col, seq_along(filled_col)))
        col_counts_df$factor_filled_col <-
            factor(
                factor_filled_col,
                levels = factor_filled_col
            )

        ggplot(
            col_counts_df,
            aes(x = "", y = col_props, fill = factor_filled_col)
        ) +
        geom_bar(stat = "identity", colour = "white", size = 0.05) +
        scale_fill_manual(
            "legend",
            values = setNames(
                object = col_counts_df$filled_col,
                as.character(col_counts_df$factor_filled_col)
            )
        ) +
        coord_polar("y", start = 0) +
        theme_void() +
        theme(legend.position = "none")
    }
}


# Set parameters ===============================================================
# Get input parameters from command line
input_integrated_rds <- snakemake@input$input_integrated_rds
ref_sample <- snakemake@params$ref_sample
output_dir <- snakemake@output$output_dir

# input_integrated_rds <- "results/Integrated_rds/Species4_Sample5.rds"
# ref_sample <- "TenX_Ptr"
# output_dir <- "results/Integrated_figures/overlap_heatmap/Species4_Sample5"


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
all_sample <- unique(seurat_umap$Sample)

# Construct MST and find subgraph centers
subtree_center <- get_mst_subtree_center(
    projection = seurat_umap[, c("UMAP_1", "UMAP_2")],
    ref_sample = ref_sample,
    all_sample = all_sample
)

# Plot contour plot of subgraph center density estimate
plotting_function_umap <- function(...) {
    center_contour_plot(
        subgraph_center = subtree_center,
        projection_umap = seurat_umap[, c("UMAP_1", "UMAP_2")],
        ref_sample = ref_sample,
        lower_percentile = 5,
        higher_percentile = 40,
        main = ""
    )
}

col_counts_df <- center_contour_plot(
    subgraph_center = subtree_center,
    projection_umap = seurat_umap[, c("UMAP_1", "UMAP_2")],
    ref_sample = ref_sample,
    lower_percentile = 5,
    higher_percentile = 40,
    main = ""
)
plotting_function_dist <- function(
    col_counts_df, type = "hist", df_path = NULL, ...
) {
    plot_col_counts(col_counts_df, type, df_path)
}

# Output figure
outname <- input_integrated_rds %>%
    strsplit(split = "/") %>%
    extract2(1) %>%
    tail(1) %>%
    sub(pattern = ".rds", replacement = "")
output_png_figure(
    plotting_function_umap,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/Overlap_", outname, "_on_", ref_sample, "_Clear.png"
        ),
    output_clear = TRUE
)
output_png_figure(
    plotting_function_umap,
    output_figure = TRUE,
    output_path =
        paste0(output_dir, "/Overlap_", outname, "_on_", ref_sample, ".png"),
    output_clear = FALSE
)

output_png_figure(
    plotting_function_dist,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/Overlap_heathist_", outname, "_on_", ref_sample, ".png"
        ),
    output_clear = FALSE,
    col_counts_df = col_counts_df,
    type = "hist"
)

output_png_figure(
    plotting_function_dist,
    output_figure = TRUE,
    output_path =
        paste0(
            output_dir,
            "/Overlap_heatpie_", outname, "_on_", ref_sample, ".png"
        ),
    output_clear = FALSE,
    output_ggplot = TRUE,
    col_counts_df = col_counts_df,
    type = "pie",
    df_path =
        paste0(
            output_dir,
            "/Overlap_heatpie_", outname, "_on_", ref_sample, ".csv"
        )
)
