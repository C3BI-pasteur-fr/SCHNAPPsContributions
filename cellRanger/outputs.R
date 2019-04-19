source("reactives.R")

myZippedReportFiles <- c("crSignificantGenes.csv")

callModule(
  pHeatMapModule,
  "crHeatmapPlotModule",
  crHeatImage
)

# output$crHeat_plot1 <- renderImage({
#   crHeatImage()
# }, deleteFile = TRUE)

# output$crPrioGenes <- DT::renderDataTable({
#   if(DEBUG)cat(file=stderr(), "output$crPrioGenes\n")
#   prioritized_genes = prioritized_genes()
#   cl5 = input$cluster5
#   if(is.null(prioritized_genes) | is.null(cl5)  )
#     return(NULL)
#
#
#   if(DEBUGSAVE)
#     save(file = "~/scShinyHubDebug/crPrioGenes.RData", list = c(ls(),ls(envir = globalenv())))
#   # load(file="~/scShinyHubDebug/crPrioGenes.RData")
#
#   dt = prioritized_genes[[cl5]]
#   return(dt)
# }, options = list(scrollX = TRUE)
# )


callModule(
  tableSelectionServer,
  "crPrioGenesTableModule",
  crPrioGenesTable
)


output$crSelectedGenes <- renderText({
  if (DEBUG) cat(file = stderr(), "crSelectedGenes\n")
  featureData <- featureDataReact()
  if (is.null(featureData)) {
    return(NULL)
  }
  top.genes <- sCA_dge()
  top.genes$symbol <-
    featureData[rownames(top.genes), "symbol"]

  paste0(top.genes$symbol[input$dge_rows_selected], ",")
})


# TODO as module ?
# cell ranger output table
output$clusters5 <- renderUI({
  if (DEBUG) cat(file = stderr(), "output$clusters\n")
  projections <- projections()
  if (is.null(projections)) {
    HTML("Please load data first")
  } else {
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    noOfClusters <- levels(as.factor(projections$dbCluster))
    selectInput(
      "cluster5",
      label = "Cluster",
      choices = c(1:noOfClusters),
      selected = 1
    )
  }
})
