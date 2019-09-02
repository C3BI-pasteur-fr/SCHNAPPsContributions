require(ElPiGraph.R)
require(plyr)

observe({
  .schnappsEnv$dimScorpiusX <- input$dimScorpiusX
  .schnappsEnv$dimScorpiusY <- input$dimScorpiusY
  .schnappsEnv$dimScorpiusCol <- input$dimScorpiusCol
  
})
# updateScorpiusInput <- reactive({
observe({  
  tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  if (is.null(.schnappsEnv$dimScorpiusX)) {
    .schnappsEnv$dimScorpiusX = colnames(tsneData)[1]
  }
  if (is.null(.schnappsEnv$dimScorpiusY)) {
    .schnappsEnv$dimScorpiusY = colnames(tsneData)[2]
  }
  if (is.null(.schnappsEnv$dimScorpiusCol)) {
    .schnappsEnv$dimScorpiusCol = "dbCluster"
  }
  updateSelectInput(session, "dimScorpiusX",
                    choices = colnames(tsneData),
                    selected = .schnappsEnv$dimScorpiusX
  )
  
  updateSelectInput(session, "dimScorpiusY",
                    choices = colnames(tsneData),
                    selected = .schnappsEnv$dimScorpiusY
  )
  updateSelectInput(session, "dimScorpiusCol",
                    choices = colnames(tsneData),
                    selected = .schnappsEnv$dimScorpiusCol
  )
  # updateNumericInput(session, "scorpMaxGenes",
  #                   value = 500
  # )
})
# elpi observers ----
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpi\n"))
  .schnappsEnv$dimElpi <- input$dimElpi
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiX\n"))
  .schnappsEnv$dimElpiX <- input$dimElpiX
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiY\n"))
  .schnappsEnv$dimElpiY <- input$dimElpiY
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiCol\n"))
  .schnappsEnv$dimElpiCol <- input$dimElpiCol
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiSeed\n"))
  .schnappsEnv$elpiSeed <- input$elpiSeed
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: ElpiMethod\n"))
  .schnappsEnv$ElpiMethod <- input$ElpiMethod
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiNumNodes\n"))
  .schnappsEnv$elpiNumNodes <- input$elpiNumNodes
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpinReps\n"))
  .schnappsEnv$elpinReps <- input$elpinReps
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiProbPoint\n"))
  .schnappsEnv$elpiProbPoint <- input$elpiProbPoint
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiStartNode\n"))
  .schnappsEnv$elpiStartNode <- input$elpiStartNode
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiEndNode\n"))
  .schnappsEnv$elpiEndNode <- input$elpiEndNode
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_num_permutations\n"))
  .schnappsEnv$elpi_num_permutations <- input$elpi_num_permutations
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_ntree\n"))
  .schnappsEnv$elpi_ntree <- input$elpi_ntree
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_ntree_perm\n"))
  .schnappsEnv$elpi_ntree_perm <- input$elpi_ntree_perm
})
observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_nGenes\n"))
  .schnappsEnv$elpi_nGenes <- input$elpi_nGenes
})
observe({
  endpoints <- traj_endpoints()
  
  # browser()
  # Can use character(0) to remove all choices
  if (is.null(endpoints)) {
    return(NULL)
  }
  # save endpoint in global variable to make sure that we don't update unnecssarily
  if (!is.null(.schnappsEnv$elpiEndpoints) & 
      length(endpoints) == length(.schnappsEnv$elpiEndpoints) &
      all(sort(endpoints) == sort(.schnappsEnv$elpiEndpoints))
  ) return(NULL)
  .schnappsEnv$elpiEndpoints = endpoints
  # Can also set the label and select items
  updateSelectInput(session, inputId = "elpiStartNode",
                    choices = endpoints,
                    selected = endpoints[1]
  )
  updateSelectInput(session, inputId = "elpiEndNode",
                    choices = endpoints,
                    selected = endpoints[length(endpoints)]
  )
})


observe({
  projections <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "dimElpiX",
                    choices = colnames(projections),
                    selected = .schnappsEnv$dimElpiX
  )
  
  updateSelectInput(session, "dimElpiY",
                    choices = colnames(projections),
                    selected = .schnappsEnv$dimElpiY
  )
  updateSelectInput(session, "dimElpiCol",
                    choices = colnames(projections),
                    selected = .schnappsEnv$dimElpiCol
  )
  updateNumericInput(session, "elpiSeed",
                     value = .schnappsEnv$elpiSeed
  )
  updateSelectInput(session, "dimElpi",
                    selected = .schnappsEnv$dimElpi
  )
  updateSelectInput(session, "ElpiMethod",
                    selected = .schnappsEnv$ElpiMethod
  )
  updateNumericInput(session, "elpiNumNodes",
                     value = .schnappsEnv$elpiNumNodes
  )
  updateNumericInput(session, "elpinReps",
                     value = .schnappsEnv$elpinReps
  )
  updateNumericInput(session, "elpiProbPoint",
                     value = .schnappsEnv$elpiProbPoint
  )
  updateNumericInput(session, "elpi_num_permutations",
                     value = .schnappsEnv$elpi_num_permutations
  )
  updateNumericInput(session, "elpi_ntree",
                     value = .schnappsEnv$elpi_ntree
  )
  updateNumericInput(session, "elpi_ntree_perm",
                     value = .schnappsEnv$elpi_ntree_perm
  )
  updateNumericInput(session, "elpi_nGenes",
                     value = .schnappsEnv$elpi_nGenes
  )
  
  
})
# observeProj ----
# update projections
observe({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "observeProj")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "observeProj")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("observeProj", id = "observeProj", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj started.\n")
  
  startNode <- input$elpiStartNode
  endNode <- input$elpiEndNode
  elpimode <- input$ElpiMethod
  psTime = traj_getPseudotime()
  scEx_log <- scEx_log()
  isolate({
    prjs <- sessionProjections$prjs
  })
  
  if (is.null(scEx_log) || is.null(psTime) || elpimode=="computeElasticPrincipalCircle") {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/observeProj.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/observeProj.RData")
  cn = paste0("traj_", startNode, "_", endNode)
  if (cn %in% colnames(prjs)) {
    return(NULL)
  }
  # browser()
  if (ncol(prjs) > 0) {
    # make sure we are working with the correct cells. This might change when cells were removed.
    prjs = prjs[colnames(scEx_log),,drop=FALSE]
    # didn't find a way to easily overwrite columns
    
    if (cn %in% colnames(prjs)) {
      prjs[, cn] <- psTime$Pt
    } else {
      prjs <- base::cbind(prjs, psTime$Pt, deparse.level = 0)
      colnames(prjs)[ncol(prjs)] <- cn
    }
    sessionProjections$prjs <- prjs
  } else {
    
    prjs <- data.frame(cn = psTime$Pt )
    rownames(prjs) = colnames(scEx_log)
    colnames(prjs)[ncol(prjs)] = cn
    sessionProjections$prjs = prjs
  }
})



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
  dimX <- .schnappsEnv$dimScorpiusX
  dimY <- .schnappsEnv$dimScorpiusY
  dimCol <- .schnappsEnv$dimScorpiusCol
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

