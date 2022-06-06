library(Matrix)
library(dplyr)
library(magrittr)
library(easyGgplot2)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
# Output figure
output_png_figure <- function(plotting_function,
                              # x, y, col, main = "",
                              output_figure = FALSE,
                              output_path = "temp.png",
                              output_clear = FALSE,
                              output_ggplot = FALSE,
                              ...) {
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

# Add parameter `scale` into `easyGgplot2::ggplot2.violinplot()`
this.ggplot2.violinplot <-
    function(
        data, xName = NULL, yName = NULL, groupName = NULL,
        addMean = FALSE, meanPointShape = 5, meanPointSize = 4,
        meanPointColor = "black", meanPointFill = "blue",
        addDot = FALSE, dotSize = 1, dotPosition = c("center", "jitter"),
        jitter = 0.2, groupColors = NULL, brewerPalette = NULL,
        ...
    ) {
    params <- list(...)
    trim <- ifelse(!is.null(params$trim), params$trim, TRUE)
    position <- ifelse(!is.null(params$position), params$position, "dodge")
    scale <- ifelse(!is.null(params$scale), params$scale, "area")
    if (is.null(yName) & !is.numeric(data)) {
        stop("yName is missing or NULL. In this case data should be a numeric vector")
    } else if (is.numeric(data)) {
        data <- cbind(y = data, x = rep(1, length(data)))
        xName <- "x"
        yName <- "y"
    }
    if (is.null(xName)) {
        data <- cbind(data, x = rep(1, nrow(data)))
        xName <- "x"
    }
    data <- data.frame(data)
    data[, xName] <- factor(data[, xName])
    if (is.null(groupName)) {
        p <- ggplot(data = data, aes_string(x = xName, y = yName))
    } else {
        p <- ggplot(
            data = data, aes_string(
            x = xName, y = yName,
            fill = groupName
        ))
        data[, groupName] <- factor(data[, groupName])
    }
    p <- p + geom_violin(trim = trim, position = position, scale = scale)
    if (addMean) {
        p <- p + stat_summary(
            fun.y = mean, geom = "point", shape = meanPointShape,
            size = meanPointSize, colour = meanPointColor, fill = meanPointFill
        )
    }
    if (addDot) {
        pms <- .dotplot_params(...)
        if (dotPosition[1] == "center") {
            p <- p + geom_dotplot(
                binaxis = "y", stackdir = "center",
                dotsize = dotSize, width = pms$width,
                stackratio = pms$stackratio,
                colour = "#7F7F7F4F"
            )
        } else {
            p <- p + geom_jitter(
                position = position_jitter(jitter),
                cex = dotSize, shape = 16,
                colour = "#7F7F7F4F"
            )
        }
    }
    if (!is.null(groupColors)) {
        p <- p + scale_fill_manual(values = groupColors)
        p <- p + scale_colour_manual(values = groupColors)
    } else if (!is.null(brewerPalette)) {
        p <- p + scale_fill_brewer(palette = brewerPalette)
        p <- p + scale_colour_brewer(
            palette = brewerPalette,
            guide = "none"
        )
    }
    p <- ggplot2.customize(p, ...)
    p
}
environment(this.ggplot2.violinplot) <- asNamespace("easyGgplot2")
assignInNamespace(
    "ggplot2.violinplot",
    this.ggplot2.violinplot,
    ns = "easyGgplot2"
)
# The call to `environment()` assures that the function will be able to call
# other hidden functions from the package.
# The call to `assignInNamespace()` assures that other functions from the
# package will call your updated version of the function.

# Set parameters ===============================================================
# Get input parameters from snakemake
log2_norm_ortho_dirs <- snakemake@input$log2_norm_ortho_dirs
plot_ortholog_csv <- snakemake@input$plot_ortholog_csv
sample_vector <- snakemake@params$sample_list
output_dir <- snakemake@output$output_dir

# log2_norm_ortho_dirs <- "results/Ortholog_UMI_matrix/Each_log2_norm/Multi_All"
# plot_ortholog_csv <- "rawdata/20220331_orthogroup_list_revised.csv"
# sample_vector <- c("TenX_Ptr", "MARSseq_Egr", "MARSseq_Tar", "TenX_Lch")
# output_dir <- "results/Ortholog_figures/Expression_boxplot"


# Implement program ============================================================
# Create output directory
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Input log2 norm ortholog UMI matrix
ortho_UMI_matrix_list <-
    sapply(
        simplify = FALSE,
        sample_vector,
        function(sample_name) {
            readRDS(paste0(log2_norm_ortho_dirs, "/", sample_name, ".rds"))
        }
    )

# Input ortholog table for plotting
plot_ortholog <- read.csv(plot_ortholog_csv)

# Wrap function for plotting
groupCol <- sapply(
    sample_vector,
    function(sample_name) {
        sp_name <- c("Ptr", "Egr", "Tar", "Lch")
        sp_col <-
            c("#D65F4780", "#D65F4780", "#00B0F080", "#FFC00080")
        for (i in seq_along(sp_name)) {
            if (grepl(sp_name[i], sample_name)) {
                return(sp_col[i])
            }
        }
        stop(paste("Unexpected species in sample", sample_name))
    }
)
plotting_function <- function(
    ortholog_cluster_no,
    ortholog_common_name,
    ...
) {
    ortholog_cluster_id <-
        paste0("Cluster_", strsplit(ortholog_cluster_no, "_")[[1]])
    ortholog_reducezero_exp_list <-
        sapply(
            simplify = FALSE,
            ortho_UMI_matrix_list,
            function(X) {
                if (any(ortholog_cluster_id %in% rownames(X))) {
                    ortholog_cluster_id <-
                        intersect(ortholog_cluster_id, rownames(X))
                    out <- apply(X[ortholog_cluster_id, , drop = FALSE], 2, sum)
                    if (sum(out == 0) >= 2) {
                        out <- c(0, 0, out[out > 0])
                    }
                } else {
                    out <- c(NA, NA)
                }
                return(out)
            }
        )
    ortholog_reducezero_exp_df <-
        data.frame(
            sample = factor(
                rep(
                    names(ortholog_reducezero_exp_list),
                    sapply(ortholog_reducezero_exp_list, length)
                ),
                levels = sample_vector
            ),
            ortholog_UMI = unlist(ortholog_reducezero_exp_list)
        )
    # ortholog_exp_list <- sapply(
    #     simplify = FALSE,
    #     ortho_UMI_matrix_list,
    #     function(X) {
    #         if (ortholog_cluster_id %in% rownames(X)) {
    #             out <- X[ortholog_cluster_id, ]
    #         } else {
    #             out <- c(NA, NA)
    #         }
    #         return(out)
    #     }
    # )
    # ortholog_exp_df <-
    #     data.frame(
    #         sample = rep(
    #             names(ortholog_exp_list),
    #             sapply(ortholog_exp_list, length)
    #         ),
    #         ortholog_UMI = unlist(ortholog_exp_list)
    #     )
    # all.equal(
    #     unlist(ortholog_exp_list),
    #     do.call(c, ortholog_exp_list)
    # ) # TRUE
    # ortholog_exp_df$sample <-
    #     factor(
    #         ortholog_exp_df$sample,
    #         levels = sample_vector
    #     )

    if (all(is.na(ortholog_reducezero_exp_df$ortholog_UMI))) {
        ggplot() +
            ggtitle(
                paste(ortholog_common_name, "is not detected in any samples.")
            )
    } else {
        violin_plot <- this.ggplot2.violinplot(
            data = ortholog_reducezero_exp_df,
            xName = "sample",
            yName = "ortholog_UMI",
            groupName = "sample",
            groupColors = groupCol,
            xShowTitle = FALSE,
            xShowTickLabel = FALSE,
            ytitle = "log2 normalized UMI",
            backgroundColor = "white",
            removePanelGrid = TRUE,
            removePanelBorder = TRUE,
            axisLine = c(0.5, "solid", "black"),
            addDot = FALSE,
            # dotSize = 1, dotPosition = "jitter", jitter = 0.2,
            scale = "count" # area, count, width
        )
        violin_plot +
            scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
            theme(legend.key.size = unit(0.5, "cm")) +
            ggtitle(ortholog_common_name)
    }
}

# Output figure
for (i_ortholog in seq(nrow(plot_ortholog))) {
    output_png_figure(
        plotting_function,
        output_figure = TRUE,
        output_path = paste0(
            output_dir, "/Ortholog_",
            plot_ortholog$OrthoGroupNo[i_ortholog],
            "_",
            gsub("/", " ", plot_ortholog$CommonName[i_ortholog]),
            ".png"
        ),
        output_clear = FALSE,
        output_ggplot = TRUE,
        ortholog_cluster_no = plot_ortholog$OrthoGroupNo[i_ortholog],
        ortholog_common_name = plot_ortholog$CommonName[i_ortholog]
    )
}
