# devtools::install_github("vqv/ggbiplot")
library(magrittr)
library(pheatmap)
library(ggbiplot)
library(VennDiagram)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================


# Set parameters ===============================================================
input_DEGtable_csvs <- file.path()
input_wholeTable_csv <- file.path()
PtrInternode_gene_abundances_csvs <- file.path()


# Implementation ===============================================================
DEG_list <-
    sapply(
        simplify = FALSE, USE.NAMES = FALSE,
        input_DEGtable_csvs,
        function(x) read.csv(x)$X
    )
sign_list <-
    sapply(
        simplify = FALSE, USE.NAMES = FALSE,
        input_DEGtable_csvs,
        function(x) {
            foo <- read.csv(x)
            stopifnot(ncol(foo) == 44)
            out <- sign(foo[, 37])
            return(out)
        }
    )
input_wholeTable <- read.csv(input_wholeTable_csv, row.names = 1)

plotting_matrix <- as.matrix(input_wholeTable[unlist(DEG_list), 1:30])

row_anno_df <- data.frame(
    row.names = rownames(plotting_matrix),
    Sign = ifelse(unlist(sign_list) == 1, "Up", "Down")
)

png(
    file.path(),
    pointsize = 12, res = 300,
    width = 28, height = 35, units = "cm"
)
pheatmap(
    plotting_matrix,
    annotation_row = row_anno_df,
    border_color = "white",
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    gaps_col = seq(6, 30, 6),
    gaps_row = cumsum(sapply(DEG_list, length)),
    cex = 1.2, lwd = 1.2
)
dev.off()


plotting_matrix <- as.matrix(input_wholeTable[, 1:30])
all_PCA <- prcomp(t(plotting_matrix), scale = TRUE, center = TRUE, retx = TRUE)
group_vector <- substr(colnames(plotting_matrix), 1, 14)
ggbiplot::ggbiplot(
    all_PCA,
    groups = group_vector,
    # labels = colnames(plotting_matrix),
    var.axes = FALSE
) +
    geom_point(aes(colour = group_vector), size = 3)
ggplot2::ggsave(
    file.path(),
    width = 20, height = 20, units = "cm",
    dpi = 300
)


PtrInternode_gene_abundances_df_list <-
    sapply(
        simplify = FALSE,
        PtrInternode_gene_abundances_csvs,
        read.csv, row.names = 1
    )
dected_genes_list <- sapply(
    PtrInternode_gene_abundances_df_list,
    function(X) rownames(X)[rowSums(X) > 0]
)
names(dected_genes_list) <- sub(
    file.path(),
    "",
    names(dected_genes_list)
)
names(dected_genes_list) <- sub(
    "_gene_abundances.csv",
    "",
    names(dected_genes_list)
)
n_genes <- sapply(dected_genes_list, length)
c(
    n_genes,
    "Union" = length(unique(unlist(dected_genes_list)))
)

venn.diagram(
    dected_genes_list,
    filename = file.path(),
    output = TRUE, imagetype = "png"
)
