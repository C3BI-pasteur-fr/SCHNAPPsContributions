require(ElPiGraph.R)
require(plyr)

# The output type has to be in line with the tablist item. I.e. plotOutput in this case
output$scropius_trajectory_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "scropius_trajectory_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scropius_trajectory_plot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scropius_trajectory_plot")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scropius_trajectory_plot", id = "scropius_trajectory_plot", duration = NULL)
  }
  
  projections <- projections()
  upI <- updateScorpiusInput() # needed to update input
  dimX <- input$dimScorpiusX
  dimY <- input$dimScorpiusY
  dimCol <- input$dimScorpiusCol
  doCalc <- input$scorpiusCalc
  
  if (is.null(projections)) {
    return(NULL)
  }
  if (!doCalc) {
    require(ggplot2)
    p1 <- ggplot(projections, aes_string(dimX, dimY, colour = dimCol)) + geom_point()
    return(p1)
  }
  
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scropius_trajectory_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/scropius_trajectory_plot.RData")
  space <- projections[, c(dimX, dimY)]
  require(SCORPIUS)
  traj <- SCORPIUS::infer_trajectory(space)
  draw_trajectory_plot(space, progression_group = projections[, dimCol], path = traj$path)
})

callModule(tableSelectionServer, "scorpiusTableMod", scorpiusModules)
# selected clusters heatmap module

scorpiusHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusHeatmapPlotReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusHeatmapPlotReactive")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusHeatmapPlotReactive")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusHeatmapPlotReactive", id = "scorpiusHeatmapPlotReactive", duration = NULL)
  }
  
  upI <- updateScorpiusInput() # needed to update input
  projections <- projections()
  traj <- scorpiusTrajectory()
  expr_sel <- scorpiusExpSel()
  modules <- scorpiusModules()
  
  dimCol <- input$dimScorpiusCol
  doCalc <- input$scorpiusCalc
  
  
  if (!doCalc | is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(traj)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusHeatmapPlot:NULL\n"))
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scorpiusHeatmapPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/scorpiusHeatmapPlot.RData")
  
  pixelratio <- session$clientData$pixelratio
  if (is.null(pixelratio)) pixelratio <- 1
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  outfile <- paste0(tempdir(), "/heatmapScorpius", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(expr_sel, traj$time, projections[, dimCol], modules,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {retVal})  
  return(retval)
  
})

callModule(
  pHeatMapModule,
  "scorpiusHeatmapPlotModule",
  scorpiusHeatmapPlotReactive
)

output$downLoadTraj <- downloadHandler(
  if (DEBUG) cat(file = stderr(), "downLoadTraj started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "downLoadTraj")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "downLoadTraj")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("downLoadTraj", id = "downLoadTraj", duration = NULL)
  }
  
  filename = paste0("scorpiusTraj.", Sys.Date(), ".csv"),
  content = function(file) {
    if (DEBUG) cat(file = stderr(), paste("downLoadTraj: \n"))
    scTRAJ <- scorpiusTrajectory()
    if (is.null(scTRAJ)) {
      return(NULL)
    }
    write.csv(scTRAJ, file)
  }
  exportTestValues(downLoadTraj = {scTRAJ})  
  return()
  
)


# --------------------------
# Elpi Graph
# --------------------------

output$elpi_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "elpi_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpi_plot")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpi_plot")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpi_plot", id = "elpi_plot", duration = NULL)
  }
  
  tree_data <- elpiTreeData()
  cep <- elpiGraphCompute()
  PointLabel <- elpiPointLabel()
  
  if (is.null(tree_data) | is.null(cep)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    base::save(file = "~/scShinyHubDebug/elpi_plot.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/elpi_plot.RData")
  
  require(ggrepel)
  p <- PlotPG(X = tree_data, TargetPG = cep[[length(cep)]], GroupsLab = PointLabel, p.alpha = 0.9)
  p[[1]] <- p[[1]] + geom_label_repel(
    data = plyr::ddply(p[[1]]$data, ~Group, summarise, meanA = mean(PCA), meanB = mean(PCB)),
    aes(x = meanA, y = meanB, label = Group),
    vjust = 1
  )
  p
})

output$elpi_histo <- renderPlot({
  if (DEBUG) cat(file = stderr(), "elpi_histo started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpi_histo")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpi_histo")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpi_histo", id = "elpi_histo", duration = NULL)
  }
  
  PointLabel <- elpiPointLabel()
  if (is.null(PointLabel)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    base::save(file = "~/scShinyHubDebug/elpi_histo.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load("~/scShinyHubDebug/elpi_histo.RData")
  
  barplot(table(PointLabel), las = 2, ylab = "Number of points")
})
