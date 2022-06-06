library(easyGgplot2)
library(reshape2)

# Load data ---------------------------
input_filename <- "Correlation_matrix_fiber.csv"
input_filename <- "Correlation_matrix_vessel.csv"
input_filename <- "Correlation_matrix_ray.csv"
data <- read.csv(input_filename, header = T)
df <- melt(data,id.vars = c("Cluster"))

# Plot data ---------------------------
## OUTPUT = 900 x 220 px
violin_plot <- ggplot2.violinplot(data = df, 
                                  xName = "Cluster", 
                                  yName = "value",
                                  xShowTitle = FALSE, 
                                  ytitle = "Correlation",
                                  xShowTickLabel = FALSE,
                                  backgroundColor = "white",
                                  removePanelGrid = TRUE,
                                  removePanelBorder = TRUE,
                                  axisLine = c(0.5, "solid", "black"),
                                  groupName = 'Cluster',
                                  groupColors = c("#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C", 
                                                  "#9467BD", "#D62728", "#E377C2", "#9B9B9B", "#4B4B4B"),
                                  faceting = TRUE, 
                                  facetingVarNames = "variable")
violin_plot + geom_violin(scale = "width") + 
              theme(legend.key.size = unit(0.5, "cm"))
 