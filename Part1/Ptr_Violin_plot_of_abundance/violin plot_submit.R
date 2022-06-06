library(reshape2)
library(easyGgplot2)

# Load data ---------------------------
input_filename <- "ANT_log2_norm_UMI_0.1.csv"
input_filename <- "AtHB8_log2_norm_UMI.csv"
input_filename <- "EXPA1_log2_norm_UMI.csv"
input_filename <- "PEAR1_log2_norm_UMI.csv"
input_filename <- "Ray_log2_norm_UMI_5.csv"
input_filename <- "VND6_log2_norm_UMI_0.1.csv"
data <- read.csv(input_filename, header = TRUE)
df <- melt(data, id.vars = c("Cluster"))

# Plot data ---------------------------
## OUTPUT = 900 x 220 px
violin_plot <- ggplot2.violinplot(data = df, 
                                  xName = "Cluster",
                                  yName = "value",
                                  xShowTitle = FALSE, 
                                  ytitle = "log2 normalized UMI",
                                  xShowTickLabel = FALSE,
                                  backgroundColor = "white",
                                  removePanelGrid = TRUE,
                                  removePanelBorder = TRUE,
                                  axisLine = c(0.5, "solid", "black"),
                                  groupName = "Cluster",
                                  groupColors = c("#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C",
                                                  "#9467BD", "#D62728", "#E377C2", "#9B9B9B", "#4B4B4B"),
                                  faceting = TRUE,
                                  facetingVarNames = "variable")
violin_plot + geom_violin(scale = "width") + 
              theme(legend.key.size = unit(0.5, "cm")) +
              scale_x_discrete(limits = c("#01", "#02", "#03", "#04", "#05", "#06", "#07", "#08")) +
              scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
