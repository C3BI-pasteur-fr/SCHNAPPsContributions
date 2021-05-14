require(shiny)
require(shinydashboard)
require(shinyBS)
require(shinycssloaders)
require(plotly)
# menu list
# defines the main entry
menuList <- list(
  menuItem("Trajectories",
           # id="trajectoryID",
           tabName = "TrajectoryList", startExpanded = FALSE,
           menuSubItem("Scorpius", tabName = "scorpiusTab"), 
           menuSubItem("ELPIGraph", tabName = "elpiGraphTab"),
           menuSubItem("Tempora", tabName = "temporaTab")
  )
)


# tabs with the actual content
tabList <- list(
  tabItem(
    "scorpiusTab",
    box(
      title = "Scorpius trajectory inference", solidHeader = TRUE, width = 12, status = "primary",
      
      fluidRow(
        column(
          width = 4,
          selectInput("dimScorpiusX",
                      label = "Component 1",
                      choices = c("tsne1", "tsne2", "tsne3"),
                      selected = "tsne1"
          ),
          numericInput("scorpRepeat",
                       label = "number of permutations for random forrest",
                       min = 1, max = 100, step = 1,
                       value = 3
          )
        ),
        column(
          width = 4,
          selectInput("dimScorpiusY",
                      label = "Component 2",
                      choices = c("tsne1", "tsne2", "tsne3"),
                      selected = "tsne2"
          ),
          numericInput("scorpMaxGenes",
                       label = "max number of Genes",
                       min = 200, max = 20000, step = 10,
                       value = 500
          )
        ),
        column(
          width = 4,
          selectInput("dimScorpiusCol",
                      label = "Color by",
                      choices = c("sampleNames", "tsne1", "tsne2", "tsne3"),
                      selected = "sampleNames"
          ),
          fileInput("trajInputFile",
                    "Choose .csv file with trajectory informaiton",
                    accept = c(
                      ".csv",
                      "text/comma-separated-values",
                      "text/tab-separated-values",
                      "text/plain",
                      ".csv",
                      ".tsv"
                    )
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          jqui_resizable(plotOutput("scropius_trajectory_plot", height = "672px"))
        )
      ),
      br(),
      fluidRow(
        column(
          width = 12, offset = 0,
          actionButton("updatetScorpiusParameters", "apply changes", width = "80%")
        )
      ),
      br(),
      # tags$h3("Heatmap "),
      fluidRow(
        column(
          width = 12,
          pHeatMapUI("scorpiusHeatmapPlotModule") # %>% withSpinner()
          # imageOutput('scorpiusHeatmapPlotModule', height = '672px')
        )
      ),
      fluidRow(
        column(
          width = 10,
          tableSelectionUi("scorpiusTableMod")
        )
      )
    )
  ),
  tabItem(
    "elpiGraphTab",
    box(
      title = "ElpiGraph trajectory inference", solidHeader = TRUE, width = 12, status = "primary",
      # tags$h3("trajectory by ElpiGraph"),
      # fluidRow(column(12,
      #                 offset = 1,
      #                 checkboxInput("elpiCalc", "calculate", FALSE)
      # )),
      fluidRow(
        column(
          3,
          selectInput(
            "dimElpiX",
            label = "Component 1",
            choices = c("tsne1", "tsne2", "tsne3"),
            selected = "tsne1"
          )
        ),
        column(
          3,
          selectInput(
            "dimElpiY",
            label = "Component 2",
            choices = c("tsne1", "tsne2", "tsne3"),
            selected = "tsne2"
          )
        ),
        column(
          3,
          selectInput(
            "dimElpiCol",
            label = "Color by",
            choices = c("dbCluster", "tsne1", "tsne2", "tsne3"),
            selected = "dbCluster"
          )
        ),
        column(
          3,
          numericInput(
            inputId = "elpiSeed",
            label = "Seed",
            value = 9,
            min = 1, max = 1000, step = 1
          )
        )
      ),
      fluidRow(
        column(
          3,
          selectInput(
            "dimElpi",
            label = "Dimensions to use",
            choices = c("elpiPCA", "components"),
            selected = "components"
          )
        ),
        column(
          3,
          selectInput(
            "ElpiMethod",
            label = "Method to use",
            choices = c(
              "computeElasticPrincipalCurve",
              "computeElasticPrincipalTree",
              "computeElasticPrincipalCircle"
            ),
            selected = "computeElasticPrincipalTree"
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpiNumNodes",
            label = "Number of nodes",
            value = 20,
            min = 10, max = 100, step = 1
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpinReps",
            label = "Number of repeats",
            value = 1,
            min = 1, max = 50, step = 1
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpiProbPoint",
            label = "probability of inclusing of a single point for each computation",
            value = 0.6,
            min = 0.1, max = 1, step = 0.1
          )
        )
      ),
      fluidRow(column(
        12,
        jqui_resizable(plotOutput("elpi_plot", height = "672px")) # %>% withSpinner()
      )),
      br(),
      fluidRow(
        column(
          width = 12, offset = 0,
          actionButton("elpiCalc", "apply changes", width = "80%")
        )
      ),
      br(),
      fluidRow(
        column(
          2,
          selectInput(
            inputId = "elpiStartNode",
            label = "start node of trajectory analysis",
            choices = c(1, 2),
            selected = 1
          )
        ),
        column(
          2,
          selectInput(
            inputId = "elpiEndNode",
            label = "end node of trajectory analysis",
            choices = c(1, 2),
            selected = 2
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpi_num_permutations",
            label = "elpi_num_permutations",
            value = 3
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpi_ntree",
            label = "elpi_ntree",
            value = 10000
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpi_ntree_perm",
            label = "elpi_ntree_perm",
            value = 1000
          )
        ),
        column(
          2,
          numericInput(
            inputId = "elpi_nGenes",
            label = "number of output genes",
            value = 50
          )
        )
      ),
      box(
        which = "plot", width = 12,
        fluidRow(column(
          12,
          # plotOutput("elpi_heatmap", height = "672px")
          pHeatMapUI("elpiHeatmapPlotModule")
        )),
        br(),
        tags$h3("table"),
        fluidRow(column(
          12,
          offset = 0,
          tableSelectionUi("elpiTableMod")
        )),
        br(),
        # tags$h3("Heatmap "),
        fluidRow(column(
          12,
          plotOutput("elpi_histo", height = "672px") # %>% withSpinner()
        ))
      )
    )
  ),
  tabItem("temporaTab",
          tabBox(title = "Tempora", width = 12, id = "temporaTab",
                 tabPanel(
                   title = "Parameters for Tempora trajectory inference", solidHeader = TRUE, width = 12, 
                   value = "temporaParameters",
                   # The id lets us use input$tabset1 on the server to find the current tab
                   id = "tabsetTempora",
                   fluidRow(column(
                     width = 12, offset = 0,
                     br("uses transformed data"))),
                   br(),
                   fluidRow(
                     column(
                       width = 12, offset = 0,
                       actionButton("updatetTemporaParameters", "apply changes", width = "80%")
                     )
                   ),
                   fluidRow(
                     column(4,
                            offset = 0,
                            selectInput("temporaCluster", "cluster points to be used",
                                        choices = c("dbCluster"),
                                        selected = defaultValue("temporaCluster", "dbCluster")
                            )
                     ),
                     column(4,
                            offset = 0,
                            selectInput("temporaFactor", "time variable",
                                        choices = c("sampleNames", "dbCluster"),
                                        selected = defaultValue("temporaFactor", "sampleNames")
                            )),
                     column(4,
                            offset = 0,
                            selectInput("temporaLevels", "Ordered time points",
                                        choices = c("AVC","MVE16.5","MV1","MV2"),
                                        multiple = TRUE,
                                        selected = defaultValue("temporaLevels", c("AVC","MVE16.5","MV1","MV2"))
                            )
                            
                            
                     )
                   ),
                   fluidRow(
                     column(
                       6,
                       offset = 0,
                       fileInput(
                         "temporaGMTFile",
                         "GMT file to use",
                         accept = c(".gmt"),
                         placeholder = "no file selected",
                         multiple = FALSE,
                       ) %>% setId(id="temporaGMTFile")
                     ),
                     column(
                       3,
                       numericInput(
                         inputId = "temporaMinSz",
                         label = "min size of gene sets",
                         value = 5
                       )
                     ),
                     column(
                       3,
                       numericInput(
                         inputId = "temporaMaxSz",
                         label = "max size of gene sets",
                         value = 200
                       )
                     )
                   ),
                   fluidRow(column(
                     4,
                     numericInput(
                       inputId = "temporaNPCs",
                       label = "N PCs to use",
                       value = 12
                     )
                   ),
                   column(
                     4,
                     numericInput(
                       inputId = "temporaDiff_thresh",
                       label = "Percent of permissible difference",
                       value = 0.01
                     )
                   ),
                   column(
                     4,
                     numericInput(
                       inputId = "temporaPval_thresh",
                       label = "max p-value",
                       value = 0.5
                     )
                   )
                   ),
                   fluidRow(column(
                     12,
                     jqui_resizable(plotOutput("tempora_screeplot", height = "672px")) # %>% withSpinner()
                   )),
                   br(),
                   actionButton("save2Hist_tempora_screeplot", "save to history"),
                   br(),
                   fluidRow(column(
                     12,
                     jqui_resizable(plotOutput("tempora_plot", height = "672px")) # %>% withSpinner()
                   )),
                   br(),
                   actionButton("save2Hist_tempora_plot", "save to history"),
                   br(),
                   checkbsTT("temporaCluster"),
                   checkbsTT("temporaFactor"),
                   checkbsTT("temporaLevels"), 
                   checkbsTT("temporaGMTFile"),
                   checkbsTT("temporaMinSz"),
                   checkbsTT("temporaMaxSz"),
                   checkbsTT("temporaNPCs"),
                   checkbsTT("temporaDiff_thresh"),
                   checkbsTT("temporaPval_thresh"),
                   # n_pcs = 12
                   # difference_threshold = 0.01
                   # pval_threshold = 0.5
                   
                 ),
                 tabPanel(
                   title = "p-Values of GO terms", solidHeader = TRUE, width = 12, 
                   value = "temporaTable",
                   id = "tabsetTemporaGOTable",
                   fluidRow(column(
                     12,
                     offset = 0,
                     tableSelectionUi("temporaGOTableMod")
                   )),
                   br(),
                   fluidRow(column(
                     12,
                     jqui_resizable(plotOutput("temporaSelectedGOs", height = "672px")) # %>% withSpinner()
                   )),
                   br(),
                   fluidRow(column(12,
                                   verbatimTextOutput("coE_temporaPWgenes"))),
                   actionButton("save2Hist_temporaSelectedGOs", "save to history")
                 ),
                 tabPanel(
                   title = "Trajectory inquiry", solidHeader = TRUE, width = 12, 
                   value = "temporaScorpius",
                   id = "tabsetTemporaScorpius",
                   fluidRow(column(
                     width = 4,
                     selectInput("dimTemporaX",
                                 label = "Component 1",
                                 choices = c("tsne1", "tsne2", "tsne3"),
                                 selected = "tsne1"
                     )),
                     column(
                       width = 4,
                       selectInput("dimTemporaY",
                                   label = "Component 2",
                                   choices = c("tsne1", "tsne2", "tsne3"),
                                   selected = "tsne1"
                       ))
                   ),
                   br(),
                   fluidRow(column(
                     12,
                     offset = 0,
                      jqui_resizable(plotly::plotlyOutput("tempora2dPlot", height = "672px")) # %>% withSpinner()
                   )),
                   br(),
                   actionButton("save2Hist_tempora2dPlot", "save to history"),
                   ######
                   br(),
                   # fluidRow(
                   #   column(
                   #     width = 12, offset = 1,
                   #     actionButton("temporaCalc", "apply changes", width = "80%")
                   #   )
                   # ),
                   # br(),
                   # fluidRow(
                   #   column(
                   #     2,
                   #     selectInput(
                   #       inputId = "temporaStartNode",
                   #       label = "start node of trajectory analysis",
                   #       choices = c(1, 2),
                   #       selected = 1
                   #     )
                   #   ),
                   #   column(
                   #     2,
                   #     selectInput(
                   #       inputId = "temporaEndNode",
                   #       label = "end node of trajectory analysis",
                   #       choices = c(1, 2),
                   #       selected = 2
                   #     )
                   #   ),
                   #   column(
                   #     2,
                   #     numericInput(
                   #       inputId = "tempora_num_permutations",
                   #       label = "tempora_num_permutations",
                   #       value = 3
                   #     )
                   #   ),
                   #   column(
                   #     2,
                   #     numericInput(
                   #       inputId = "tempora_ntree",
                   #       label = "tempora_ntree",
                   #       value = 10000
                   #     )
                   #   ),
                   #   column(
                   #     2,
                   #     numericInput(
                   #       inputId = "tempora_ntree_perm",
                   #       label = "tempora_ntree_perm",
                   #       value = 1000
                   #     )
                   #   ),
                   #   column(
                   #     2,
                   #     numericInput(
                   #       inputId = "tempora_nGenes",
                   #       label = "number of output genes",
                   #       value = 50
                   #     )
                   #   )
                   # ),
                   # box(
                   #   which = "plotTempora", width = 12,
                   #   fluidRow(column(
                   #     12,
                   #     # plotOutput("tempora_heatmap", height = "672px")
                   #     pHeatMapUI("temporaHeatmapPlotModule")
                   #   )),
                   #   br(),
                   #   tags$h3("table"),
                   #   fluidRow(column(
                   #     10,
                   #     offset = 1,
                   #     tableSelectionUi("temporaTableMod")
                   #   )),
                   #   br()
                   # )
                 )
                 
                 
          )
  )
)
# declare heavy calculations
# myHeavyCalculations <- NULL
