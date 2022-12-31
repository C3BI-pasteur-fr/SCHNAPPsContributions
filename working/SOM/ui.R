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
                   sc_numericInput("coE_dimSOM", "number of nodes per dimension",
                                   defaultValue("coE_dimSOM", 20),
                                min = 2, max = 100
                   )
            ), 
            column(width = 3,
                   sc_textInput("coE_geneSOM", "Gene of interest", value = defaultValue("coE_geneSOM", "notyet"))
            ),
            column(width = 3,
                   sc_selectInput("coE_distSOM",
                               label = "Distance",
                               choices = c("raw", "Spearman", "standardized"),
                               selected = defaultValue("coE_distSOM","raw")
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
        sc_numericInput("coE_dimSOMX", "row",
                     defaultValue("coE_dimSOMX",1),
                     min = 1, max = 100
        )
      ),
      column(
        width = 6,
        sc_numericInput("coE_dimSOMY", "column",
                     defaultValue("coE_dimSOMY",1),
                     min = 1, max = 100
        )
      ),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somInfo")
      )),
      br(),
      fluidRow(column(
        width = 12,
        verbatimTextOutput("coE_somInfoSymbol")
      ))
      )
    )
  )
)