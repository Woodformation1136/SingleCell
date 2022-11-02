library(shiny)

# Define UI ----
ui <- fluidPage(

    # App title ----
    titlePanel("Plot transcript abundance on UMAP"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(

            h4(strong("1. Set the plotting data")),

            radioButtons(
                inputId = "UMAP_type",
                label = "Choose the source of UMAP:",
                choices = c(
                    "single sample analysis (Cellranger)",
                    "multiple samples analysis (Seurat CCA)"
                ),
                selected = "multiple samples analysis (Seurat CCA)"
            ),

            uiOutput("multisample_selection"),

            radioButtons(
                inputId = "UMI_type",
                label = "Choose the type of UMI counts:",
                choices = c(
                    "gene UMI counts",
                    "ortholog UMI counts"
                ),
                selected = "ortholog UMI counts"
            ),

            uiOutput("plotting_sample_selection")
        ),

        # Main panel for displaying outputs ----
        mainPanel(

            h4(strong("2. Collect plotting data (UMI data)")),

            strong("UMI data in used:"),

            actionButton(
                inputId = "load_UMI_data",
                label = "Update"
            ),

            verbatimTextOutput("UMI_data_in_use", placeholder = TRUE),


            h4(strong("3. Set the plotting parameters")),

            uiOutput("plotting_feature_selection"),
            uiOutput("plotting_feature_selection_eg"),

            p(strong("Check to output clear figure:")),
            checkboxInput(
                inputId = "output_clear",
                label = "",
                value = FALSE
            ),

            sliderInput(
                inputId = "color_range",
                label = "Set the color gradient range:",
                min = 0,
                max = 100,
                value = c(20, 50),
                step = 1,
                round = TRUE,
                ticks = FALSE
            ),

            actionButton(
                inputId = "figure_out",
                label = "Figure out"
            ),

            imageOutput(
                outputId = "image"
            ),


            h4(strong("4. Download the figure")),

            uiOutput("download_filename_setting"),

            downloadButton("download_figure", "Download")
        )
    )
)