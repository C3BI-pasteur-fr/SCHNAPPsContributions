require(shiny)
require(shinydashboard)
require(shinyBS)
require(shinycssloaders)
# menu list
# defines the main entry
menuList <- list(
  menuItem("Trajectories",
           # id="trajectoryID",
    tabName = "TrajectoryList", startExpanded = FALSE,
    menuSubItem("Scorpius", tabName = "scorpiusTab"), menuSubItem("ELPIGraph", tabName = "elpiGraphTab")
  )
)


# tabs with the actual content
tabList <- list(
  crHeatMapTab =
    tabItem(
      "scorpiusTab",
      tags$h3("trajectory in 2D space"),
      fluidRow(column(12,
        offset = 1,
        checkboxInput("scorpiusCalc", "calculate", FALSE)
      )),
      fluidRow(
        column(
          3,
          selectInput(
            "dimScorpiusX",
            label = "Component 1",
            choices = c("tsne1", "tsne2", "tsne3"),
            selected = "tsne1"
          )
        ),
        column(
          3,
          selectInput(
            "dimScorpiusY",
            label = "Component 2",
            choices = c("tsne1", "tsne2", "tsne3"),
            selected = "tsne2"
          )
        ),
        column(
          3,
          selectInput(
            "dimScorpiusCol",
            label = "Color by",
            choices = c("sample", "tsne1", "tsne2", "tsne3"),
            selected = "sample"
          )
        ),
        column(
          3,
          numericInput(
            "scorpMaxGenes",
            label = "max number of Genes",
            min = 200, max = 20000, step = 10,
            value = 500
          )
        )
      ),
      fluidRow(      
        column(
        3,
        numericInput(
          "scorpRepeat",
          label = "number of permutations for random forrest",
          min = 1, max = 100, step = 1,
          value = 3
        )
      )
      ),
      fluidRow(column(
        12,
        tipify(
          downloadButton("downLoadTraj", "Download trajectory"),
          "<h3>download trajectory to csv file</h3>"
        ), fileInput(
          "trajInputFile",
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
      )),
      fluidRow(column(
        12,
        plotOutput("scropius_trajectory_plot", height = "672px") # %>% withSpinner()
      )),
      # tags$h3("Heatmap "),
      fluidRow(column(
        12,
        pHeatMapUI("scorpiusHeatmapPlotModule") %>% withSpinner()
        # imageOutput('scorpiusHeatmapPlotModule', height = '672px') #%>% withSpinner()
      )),
      tags$h3("table"),
      fluidRow(column(
        10,
        offset = 1,
        tableSelectionUi("scorpiusTableMod")
      ))
    ),
  elpiTab = tabItem(
    "elpiGraphTab",
    tags$h3("trajectory by ElpiGraph"),
    fluidRow(column(12,
      offset = 1,
      checkboxInput("elpiCalc", "calculate", FALSE)
    )),
    fluidRow(
      column(
        4,
        selectInput(
          "dimElpi",
          label = "Dimensions to use",
          choices = c("elpiPCA"),
          selected = "elpiPCA"
        )
      ),
      column(
        4,
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
        4,
        numericInput(
          inputId = "elpiNumNodes",
          label = "Number of nodes",
          value = 20,
          min = 10, max = 100, step = 1
        )
      ),
      column(
        4,
        numericInput(
          inputId = "elpinReps",
          label = "Number of repeats",
          value = 1,
          min = 1, max = 50, step = 1
        )
      ),
      column(
        4,
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
      plotOutput("elpi_plot", height = "672px") # %>% withSpinner()
    )),
    fluidRow(column(
      4,
      numericInput(
        inputId = "elpiStartNode",
        label = "start node of trajectory analysis",
        value = 1
      ),
      numericInput(
        inputId = "elpiEndNode",
        label = "end node of trajectory analysis",
        value = 2
      ),
      numericInput(
        inputId = "elpi_num_permutations",
        label = "elpi_num_permutations",
        value = 3
      ),
      numericInput(
        inputId = "elpi_ntree",
        label = "elpi_ntree",
        value = 10000
      ),
      numericInput(
        inputId = "elpi_ntree_perm",
        label = "elpi_ntree_perm",
        value = 1000
      ),
      numericInput(
        inputId = "elpi_nGenes",
        label = "number of output genes",
        value = 50
      )
    )),
    fluidRow(column(
      12,
      plotOutput("elpi_heatmap", height = "672px")
    )),
    # tags$h3("Heatmap "),
    fluidRow(column(
      12,
      plotOutput("elpi_histo", height = "672px") # %>% withSpinner()
    ))
  )
)

# declare heavy calculations
myHeavyCalculations <- NULL
