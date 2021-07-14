library(dplyr)
library(reshape2)
library(easyGgplot2)

# Load data ---------------------------
input_filename <- "Cluster_cellnumber.csv"
data <- read.csv(input_filename, header = TRUE)
df <- melt(data, id.vars = c("Cluster"))

Ptr <- filter(df, variable == "Ptr")
Egr <- filter(df, variable == "Egr")
Lch <- filter(df, variable == "Lch")
Tar <- filter(df, variable == "Tar")

# Plot data ---------------------------
## OUTPUT = 720 X 200px
### Ptr
Ptr_plot <- ggplot2.barplot(data = Ptr,
                            xName = "Cluster", 
                            yName = "value",
                            xShowTitle = FALSE, 
                            yShowTitle = FALSE,
                            orientation = "horizontal",
                            xShowTickLabel = TRUE,
                            backgroundColor = "white",
                            removePanelGrid = TRUE,
                            removePanelBorder = TRUE,
                            groupName = "Cluster",
                            groupColors = c("#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C",
                                            "#9467BD", "#D62728", "#E377C2", "#9B9B9B", "#4B4B4B"),
                            legendItemOrder = c("#10", "#09", "#08", "#07", "#06",
                                                "#05", "#04", "#03", "#02", "#01"),
                            axisLine = c(0.5, "solid", "black"))
Ptr_plot + geom_text(aes(label = value), vjust = 0.4, hjust = -0.3) + 
           theme(legend.position = "none") +
           scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Egr
Egr_plot <- ggplot2.barplot(data = Egr,
                            xName = "Cluster", 
                            yName = "value",
                            xShowTitle = FALSE, 
                            yShowTitle = FALSE,
                            orientation = "horizontal",
                            xShowTickLabel = TRUE,
                            backgroundColor = "white",
                            removePanelGrid = TRUE,
                            removePanelBorder = TRUE,
                            groupName = "Cluster",
                            groupColors = c("#9467BD", "#F8E71C", "#8C564B", "#1F77B4", "#2AA02A",
                                            "#FF7F0F", "#D62728", "#E377C2", "#4B4B4B", "#9B9B9B"),
                            legendItemOrder = c("#09", "#10", "#08", "#07", "#01",
                                                "#02", "#05", "#06", "#03", "#04"),
                            axisLine = c(0.5, "solid", "black"))
Egr_plot + geom_text(aes(label = value), vjust = 0.4, hjust = -0.3) + 
           theme(legend.position = "none") +
           scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Tar
Tar_plot <- ggplot2.barplot(data = Tar, 
                            xName = "Cluster", 
                            yName = "value",
                            xShowTitle = FALSE, 
                            yShowTitle = FALSE,
                            orientation = "horizontal",
                            xShowTickLabel = TRUE,
                            backgroundColor = "white",
                            removePanelGrid = TRUE,
                            removePanelBorder = TRUE,
                            groupName = "Cluster", 
                            groupColors = c("#1F77B4", "#E377C2", "#55A3FF", "#8C564B", "#9467BD",
                                          "#FF7F0F", "#46CBE5", "#F8E71C", "#2AA02A", "#9B9B9B"),
                            legendItemOrder = c("#10", "#03", "#07", "#02", "#05",
                                                "#08", "#09", "#06", "#04", "#01"),
                            axisLine = c(0.5, "solid", "black"))
Tar_plot + geom_text(aes(label = value), vjust = 0.4, hjust = -0.3) +
           theme(legend.position = "none") +
           scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Lch
Lch_plot <- ggplot2.barplot(data = Lch, 
                            xName = "Cluster", 
                            yName = "value",
                            xShowTitle = FALSE, 
                            yShowTitle = FALSE,
                            orientation = "horizontal",
                            xShowTickLabel = TRUE,
                            backgroundColor = "white",
                            removePanelGrid = TRUE,
                            removePanelBorder = TRUE,
                            groupName = "Cluster",
                            groupColors = c("#1F77B4", "#8C564B", "#d62728", "#ff7f0f", "#f8e71c",
                                            "#e377c2", "#9467bd", "#2aa02a", "#63ee9b", "#9B9B9B"),
                            legendItemOrder = c("#10", "#09", "#06", "#03", "#07", 
                                                "#05", "#08", "#04", "#02", "#01"),
                            axisLine = c(0.5, "solid", "black"))
Lch_plot + geom_text(aes(label = value), vjust = 0.4, hjust = -0.3) +
           theme(legend.position = "none") +
           scale_y_continuous(expand = expansion(mult = c(0, .1)))
