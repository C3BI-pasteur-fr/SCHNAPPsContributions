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
  
  traj <- scorpiusTrajectory()
  projections <- projections()
  space <- scorpiusSpace()
  # upI <- updateScorpiusInput() # needed to update input
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
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scropius_trajectory_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scropius_trajectory_plot.RData")
  # space <- projections[, c(dimX, dimY)]
  require(SCORPIUS)
  # traj <- SCORPIUS::infer_trajectory(space)
  colnames(traj) = c("Comp1", "Comp2", "time")
  draw_trajectory_plot(space, progression_group = projections[rownames(space), dimCol], path = as.matrix(traj[,1:2]))
})

callModule(tableSelectionServer, "scorpiusTableMod", scorpiusModulesTable)
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
  
  # upI <- updateScorpiusInput() # needed to update input
  projections <- projections()
  traj <- scorpiusTrajectory()
  expr_sel <- scorpiusExpSel()
  modules <- scorpiusModules()
  
  dimCol <- input$dimScorpiusCol
  doCalc <- input$scorpiusCalc
  pixelratio <- session$clientData$pixelratio
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  if (!doCalc | is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(traj)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusHeatmapPlot:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusHeatmapPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusHeatmapPlot.RData")
  
  if (is.null(pixelratio)) pixelratio <- 1
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  outfile <- paste0(tempdir(), "/heatmapScorpius", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(expr_sel$expr_sel, traj$time, projections[rownames(expr_sel$expr_sel), dimCol], modules,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {retVal})  
  return(retVal)
  
})

callModule(
  pHeatMapModule,
  "scorpiusHeatmapPlotModule",
  scorpiusHeatmapPlotReactive
)

output$downLoadTraj <- downloadHandler(
  filename = paste0("scorpiusTraj.", Sys.Date(), ".csv"),
  content = function(file) {
    if (DEBUG) cat(file = stderr(), paste("downLoadTraj: \n"))
    traj <- scorpiusTrajectory()
    if (is.null(traj)) {
      return(NULL)
    }
    if (.schnappsEnv$DEBUGSAVE) {
      save(file = "~/SCHNAPPsDebug/downLoadTraj.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/SCHNAPPsDebug/downLoadTraj.RData")
    write.csv(traj, file)
  }
)


# --------------------------
# Elpi Graph
# --------------------------


elpiHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "elpiHeatmapPlotReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiHeatmapPlotReactive")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpiHeatmapPlotReactive")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpiHeatmapPlotReactive", id = "elpiHeatmapPlotReactive", duration = NULL)
  }
  
  # upI <- updateScorpiusInput() # needed to update input
  projections <- projections()
  psTime = traj_getPseudotime()
  expr_sel <- traj_elpi_gimp()
  modules <- traj_elpi_modules()
  
  dimCol <- input$dimElpiCol
  doCalc <- input$elpiCalc
  pixelratio <- session$clientData$pixelratio
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/elpiHeatmapPlotReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/elpiHeatmapPlotReactive.RData")

    if (!doCalc | is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(psTime)) {
    if (.schnappsEnv$DEBUG) cat(file = stderr(), paste("scorpiusHeatmapPlot:NULL\n"))
    return(NULL)
  }
  
  if (is.null(pixelratio)) pixelratio <- 1
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  outfile <- paste0(tempdir(), "/heatmapScorpius", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  
  pst = psTime$Pt[which(!is.na(psTime$Pt))]
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(expr_sel$expr_sel, time = pst,  projections[rownames(expr_sel$expr_sel), dimCol], modules,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {retVal})  
  return(retVal)
  
})

callModule(
  pHeatMapModule,
  "elpiHeatmapPlotModule",
  elpiHeatmapPlotReactive
)

# 
# output$elpi_heatmap <- renderPlot({
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "traj_getPseudotime")
#     if (!is.null(getDefaultReactiveDomain()))
#       removeNotification(id = "traj_getPseudotime")
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("traj_getPseudotime", id = "traj_getPseudotime", duration = NULL)
#   }
#   if (.schnappsEnv$DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
#   scEx_log <- scEx_log()
#   projections <- projections()
#   TreeEPG <- elpiGraphCompute()
#   elpimode <- input$ElpiMethod
#   tree_data <- elpiTreeData()
#   tragetPath <- traj_tragetPath()
#   gene_sel <- traj_elpi_gimp() 
#   modules <-  traj_elpi_modules()
#   psTime = traj_getPseudotime()
#   
#   if (is.null(gene_sel) || is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle" || is.null(psTime)) {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/traj_getPseudotime.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/SCHNAPPsDebug/traj_getPseudotime.RData")
#   
#   ## Select most important genes (set ntree to at least 10000!)
#   # gene_sel <- geneImport[1:50,]
#   gene_sel = gene_sel$gene_sel
#   expr_sel <- t(as.matrix(assays(scEx_log)[[1]][gene_sel$gene,which(!is.na(psTime$Pt))]))
#   
#   pst = psTime$Pt[which(!is.na(psTime$Pt))]
#   
#   
#   p <- SCORPIUS::draw_trajectory_heatmap(x = expr_sel, time = pst, progression_group = projections$dbCluster[which(!is.na(psTime$Pt))] ,
#                                          modules=modules, show_labels_row = TRUE)
#   
# })

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
  doCalc <- input$elpiCalc
  projections <- projections()
  dimX <- input$dimElpiX
  dimY <- input$dimElpiY
  dimCol <- input$dimElpiCol
  
  if (!doCalc) {
    require(ggplot2)
    p1 <- ggplot(projections, aes_string(dimX, dimY, colour = dimCol)) + geom_point()
    return(p1)
  }
  
  if (is.null(tree_data) | is.null(cep)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpi_plot.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/elpi_plot.RData")
  
  require(ggrepel)
  
  NodeLabs <- 1:nrow(cep[[length(cep)]]$NodePositions)
  NodeLabs[degree(ConstructGraph(cep[[length(cep)]])) != 1] <- NA
  
  p <- PlotPG(X = tree_data, TargetPG = cep[[length(cep)]],
              NodeLabels = NodeLabs,
              LabMult = 5, PointSize = NA, p.alpha = .1)
  
  
  # p <- PlotPG(X = tree_data, TargetPG = cep[[length(cep)]], GroupsLab = PointLabel, p.alpha = 0.9)
  # p[[1]] <- p[[1]] + geom_label_repel(
  #   data = plyr::ddply(p[[1]]$data, ~Group, summarise, meanA = mean(PCA), meanB = mean(PCB)),
  #   aes(x = meanA, y = meanB, label = Group),
  #   vjust = 1
  # )
  p
})

output$elpi_moduleHeatmap <- renderPlot({
  draw_trajectory_heatmap(expr_sel, pst, group_name, modules)
})

callModule(tableSelectionServer, "elpiTableMod", elpiModulesTable)

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
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpi_histo.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load("~/SCHNAPPsDebug/elpi_histo.RData")
  
  barplot(table(PointLabel), las = 2, ylab = "Number of points")
})

