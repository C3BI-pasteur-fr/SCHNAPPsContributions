library(CelliD)
menuList <- list(
  shinydashboard::menuItem("Cell-ID",
                           icon = icon("dashboard"),
                           # id="coexpressionID",
                           tabName = "cellID", startExpanded = FALSE,
                           shinydashboard::menuSubItem("cellID", tabName = "cellIDTab")
  )
)

# cellIDTab ----
tabList = list(
  shinydashboard::tabItem(
    "cellIDTab",
    tabBox(title = "Cell-ID", width = 12, id = "cellIDMain",
           tabPanel(
             title = "CellID params", width = 12, value = "cellIDParam",
             # footer = "see https://github.com/RausellLab/CelliD for further information",
             fluidRow(
               column(
                 width = 12, offset = 1,
                 actionButton("updateCellID", "apply changes", width = "80%")
               )
             ),
             br(),
             fluidRow(
               # cellSelectionUI("cellID_Method"),
               box(width = 12,
                   fluidRow(
                     column(width = 3,
                            sc_selectizeInput("cellID_Method", "Method to use",
                                           choices = c("MCA", "DMAP", "TSNE", "UMAP"), selected = defaultValue("cellID_Method", 'MCA')
                            )
                     ), 
                     column(width = 3, 
                            sc_numericInput(inputId = "cellID_nmcs", label = "number of components", 
                                            value = defaultValue("cellID_nmcs",20), min = 2)), 
                     column(width = 3, 
                            sc_selectizeInput(inputId = "cellID_gBy", label = "factor to use for grouping", 
                                           choices = c("dbCluster"), selected = defaultValue("cellID_gBy", "dbCluster")))
                   )
               ),
               box(width = 12,
                   fluidRow(
                     column(width = 12,
                            tableSelectionUi("cellID_GeneSetTable")
                     ))
               )
             )
           ),
           
           tabPanel(title = "cellTable", width = 12, id = "cellIDcellTable",
                    fluidRow(
                      column(
                        width = 12, offset = 0,
                        tableSelectionUi("cellID_cellTable")
                      )
                    )
           ),
           tabPanel(title = "geneTable", width = 12, id = "cellIDgeneTable",
                    fluidRow(
                      column(
                        width = 12, offset = 0,
                        tableSelectionUi("cellID_geneTable")
                      )
                    )
           ),
           tabPanel(title = "stdevTable", width = 12, id = "cellIDstdevTable",
                    fluidRow(
                      column(
                        width = 12, offset = 0,
                        # tableSelectionUi("cellID_stdevTable")
                        plotOutput("cellID_stdev")
                      )
                    )
           )
    )
  )
  
)