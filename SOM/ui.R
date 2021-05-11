menuList <- list(
  shinydashboard::menuItem("SOM",
                           icon = icon("dashboard"),
                           # id="coexpressionID",
                           tabName = "som", startExpanded = FALSE,
                           shinydashboard::menuSubItem("SOM cluster", tabName = "SOMcluster")
  )
)

# SOMcluster ----
tabList = list(
  shinydashboard::tabItem(
    "SOMcluster",
    box(
      title = "Self organizing map (SOM)", solidHeader = TRUE, width = 12, status = "primary",
      footer = "Here, we calculate a SOM on all genes using the information from all cells. Then we ask, which other genes are in the same cluster as the gene of intereset.",
      fluidRow(
        column(
          width = 12, offset = 1,
          actionButton("updateSOMParameters", "apply changes", width = "80%")
        )
      ),
      br(),
      fluidRow(
        cellSelectionUI("coE_SOM_dataInput"),
        box(
          fluidRow(
            column(width = 3,
                   numericInput("coE_dimSOM", "number of nodes per dimension",
                                20,
                                min = 2, max = 100
                   )
            ), 
            column(width = 3,
                   textInput("coE_geneSOM", "Gene of interest", value = defaultValueSingleGene)
            ),
            column(width = 3,
                   selectInput("coE_distSOM",
                               label = "Distance",
                               choices = c("raw", "Spearman", "standardized"),
                               selected = "raw"
                   )
            )
          )
        ),
        
        # column(width = 3,
        #        selectInput(inputId = "coE_clusterSOM", label = "Clusters/Factor to use", 
        #                    choices = c("dbCluster", "sampleName"),
        #                    selected = "dbCluster")
        # ),
        # column(width = 3,
        #        selectInput(inputId = "coE_clusterValSOM", label = "Values to use",
        #                    choices = c("1","2"), selected = "1", multiple = TRUE)
        # )
      ),
      br(),
      fluidRow(column(
          width = 12,
          pHeatMapUI("coE_heatmapSOM")
        )
      ),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somGenes")
      )),
      br(),
      fluidRow(column(
        width = 12,
        plotOutput("coE_SOMcodebook")
      )),
      br(),
      fluidRow(column(
        width = 12,
        plotOutput("coE_SOMcomponents")
      )),
      br(),
      fluidRow(column(
        width = 12,
        plotOutput("coE_SOMuMat")
      )),
      br(),
      fluidRow(column(
        width = 6,
        numericInput("coE_dimSOMX", "row",
                     1,
                     min = 1, max = 100
        )
      ),
      column(
        width = 6,
        numericInput("coE_dimSOMY", "column",
                     1,
                     min = 1, max = 100
        )
      ),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somInfo")
      ))
      )
    )
  )
)