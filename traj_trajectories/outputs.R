# Loading required package: plyr
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   You have loaded plyr after dplyr - this is likely to cause problems.
# If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
#   library(plyr); library(dplyr)
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
require(ElPiGraph.R)
require(plyr)
library(dplyr)



# observe scorpius Proj ----
# update projections: add new projections
observe(label = "ob20sc", {
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "observeScorpiusProj")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "observeScorpiusProj")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("observeScorpiusProj", id = "observeScorpiusProj", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeScorpiusProj 20sc started.\n")
  
  traj <- scorpiusTrajectory()
  scEx_log <- scEx_log()
  isolate({
    prjs <- sessionProjections$prjs
  })
  
  if (is.null(scEx_log) || is.null(traj)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/observeScorpiusProj.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/observeScorpiusProj.RData")
  cn <- "traj_scorpius"
  # if (cn %in% colnames(prjs)) {
  #   return(NULL)
  # }
  # browser()
  if (ncol(prjs) > 0) {
    # make sure we are working with the correct cells. This might change when cells were removed.
    prjs <- prjs[colnames(scEx_log), , drop = FALSE]
    # didn't find a way to easily overwrite columns
    
    if (cn %in% colnames(prjs)) {
      prjs[, cn] <- traj$time
    } else {
      prjs <- base::cbind(prjs, traj$time, deparse.level = 0)
      colnames(prjs)[ncol(prjs)] <- cn
    }
    sessionProjections$prjs <- prjs
  } else {
    prjs <- data.frame(cn = traj$time)
    rownames(prjs) <- colnames(scEx_log)
    colnames(prjs)[ncol(prjs)] <- cn
    sessionProjections$prjs <- prjs
  }
})

observe(label = "ob_scorpButton",
        {
          if (DEBUG) cat(file = stderr(), "observe ob_scorpButton\n")
          input$updatetScorpiusParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("scorpDimX", input$dimScorpiusX),
              c("scorpDimY", input$dimScorpiusY),
              c("scorpMaxGenes", input$scorpMaxGenes),
              c("scorpRepeat", input$scorpRepeat),
              c("scorpInFile", input$trajInputFile)
            )
          )
          
          updateButtonColor(buttonName = "updatetScorpiusParameters", parameters = c(
            "scorpDimX", "scorpDimY",
            "scorpMaxGenes", "scorpRepeat" , "scorpInFile"
          ))
        })

observe(label = "ob1", {
  .schnappsEnv$dimScorpiusX <- input$dimScorpiusX
  .schnappsEnv$dimScorpiusY <- input$dimScorpiusY
  .schnappsEnv$dimScorpiusCol <- input$dimScorpiusCol
})
# updateScorpiusInput <- reactive({
observe(label = "ob2", {
  tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  if (is.null(.schnappsEnv$dimScorpiusX)) {
    .schnappsEnv$dimScorpiusX <- colnames(tsneData)[1]
  }
  if (is.null(.schnappsEnv$dimScorpiusY)) {
    .schnappsEnv$dimScorpiusY <- colnames(tsneData)[2]
  }
  if (is.null(.schnappsEnv$dimScorpiusCol)) {
    .schnappsEnv$dimScorpiusCol <- "dbCluster"
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
observe(label = "ob3", {
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpi\n"))
  .schnappsEnv$dimElpi <- input$dimElpi
})
observe(label = "ob4", {
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiX\n"))
  .schnappsEnv$dimElpiX <- input$dimElpiX
})
observe(label = "ob5", {
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiY\n"))
  .schnappsEnv$dimElpiY <- input$dimElpiY
})
observe(label = "ob6", {
  if (DEBUG) cat(file = stderr(), paste0("observe: dimElpiCol\n"))
  .schnappsEnv$dimElpiCol <- input$dimElpiCol
})
observe(label = "ob7", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiSeed\n"))
  .schnappsEnv$elpiSeed <- input$elpiSeed
})
observe(label = "ob8", {
  if (DEBUG) cat(file = stderr(), paste0("observe: ElpiMethod\n"))
  .schnappsEnv$ElpiMethod <- input$ElpiMethod
})
observe(label = "ob9", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiNumNodes\n"))
  .schnappsEnv$elpiNumNodes <- input$elpiNumNodes
})
observe(label = "ob10", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpinReps\n"))
  .schnappsEnv$elpinReps <- input$elpinReps
})
observe(label = "ob11", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiProbPoint\n"))
  .schnappsEnv$elpiProbPoint <- input$elpiProbPoint
})
observe(label = "ob12", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiStartNode", input$elpiStartNode, "\n"))
  .schnappsEnv$elpiStartNode <- input$elpiStartNode
})
observe(label = "ob13", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpiEndNode:", input$elpiEndNode,"\n"))
  .schnappsEnv$elpiEndNode <- input$elpiEndNode
})
observe(label = "ob14", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_num_permutations\n"))
  .schnappsEnv$elpi_num_permutations <- input$elpi_num_permutations
})
observe(label = "ob15", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_ntree\n"))
  .schnappsEnv$elpi_ntree <- input$elpi_ntree
})
observe(label = "ob16", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_ntree_perm\n"))
  .schnappsEnv$elpi_ntree_perm <- input$elpi_ntree_perm
})
observe(label = "ob17", {
  if (DEBUG) cat(file = stderr(), paste0("observe: elpi_nGenes\n"))
  .schnappsEnv$elpi_nGenes <- input$elpi_nGenes
})

# observe traj_endpoints ----
observe(label = "ob18", {
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
  ) {
    return(NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj 18 started.\n")
  
  .schnappsEnv$elpiEndpoints <- endpoints
  # Can also set the label and select items
  updateSelectInput(session,
                    inputId = "elpiStartNode",
                    choices = endpoints,
                    selected = endpoints[1]
  )
  updateSelectInput(session,
                    inputId = "elpiEndNode",
                    choices = endpoints,
                    selected = endpoints[length(endpoints)]
  )
})

# observe projections ----
observe(label = "ob19", {
  projections <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj 19 started.\n")
  
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
# update projections: add new projections
observe(label = "ob20", {
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "observeProj")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "observeProj")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("observeProj", id = "observeProj", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "observeProj 20 started.\n")
  
  startNode <- input$elpiStartNode
  endNode <- input$elpiEndNode
  elpimode <- input$ElpiMethod
  psTime <- traj_getPseudotime()
  scEx_log <- scEx_log()
  isolate({
    prjs <- sessionProjections$prjs
  })
  
  if (is.null(scEx_log) || is.null(psTime) || elpimode == "computeElasticPrincipalCircle") {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/observeProj.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/observeProj.RData")
  cn <- paste0("traj_", startNode, "_", endNode)
  if (cn %in% colnames(prjs)) {
    return(NULL)
  }
  # browser()
  if (ncol(prjs) > 0) {
    # make sure we are working with the correct cells. This might change when cells were removed.
    prjs <- prjs[colnames(scEx_log), , drop = FALSE]
    # didn't find a way to easily overwrite columns
    
    if (cn %in% colnames(prjs)) {
      prjs[, cn] <- psTime$Pt
    } else {
      prjs <- base::cbind(prjs, psTime$Pt, deparse.level = 0)
      colnames(prjs)[ncol(prjs)] <- cn
    }
    sessionProjections$prjs <- prjs
  } else {
    prjs <- data.frame(cn = psTime$Pt)
    rownames(prjs) <- colnames(scEx_log)
    colnames(prjs)[ncol(prjs)] <- cn
    sessionProjections$prjs <- prjs
  }
})

# set scorpiuseParameters on button pressed----
observeEvent(input$updatetScorpiusParameters,{
  # only react on this value
  cat(file = stderr(), paste("hit button\n"))
  input$updatetScorpiusParameters
  scorpiuseParameters$dimX = isolate(input$dimScorpiusX)
  scorpiuseParameters$dimY = isolate(input$dimScorpiusY)
  scorpiuseParameters$scInput = isolate(scorpiusInput())
  scorpiuseParameters$scorpMaxGenes = isolate(input$scorpMaxGenes)
  scorpiuseParameters$scorpRepeat = isolate(input$scorpRepeat)
  
  setRedGreenButton(
    vars = list(
      c("scorpDimX", scorpiuseParameters$dimX),
      c("scorpDimY", scorpiuseParameters$dimY),
      c("scorpMaxGenes", isolate(input$scorpMaxGenes)),
      c("scorpRepeat", isolate(input$scorpRepeat)),
      c("scorpInFile", isolate(input$trajInputFile))
    ),
    button = "updatetScorpiusParameters"
  )
  
})


# scropius_trajectory_plot ----
# The output type has to be in line with the tablist item. I.e. plotOutput in this case
output$scropius_trajectory_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "scropius_trajectory_plot started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scropius_trajectory_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scropius_trajectory_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scropius_trajectory_plot", id = "scropius_trajectory_plot", duration = NULL)
  }
  
  traj <- scorpiusTrajectory()
  projections <- isolate(projections())
  space <- isolate(scorpiusSpace())
  # upI <- updateScorpiusInput() # needed to update input
  dimX <- input$dimScorpiusX
  dimY <- input$dimScorpiusY
  dimCol <- input$dimScorpiusCol
  # doCalc <- input$scorpiusCalc
  
  if (is.null(projections) ) {
    if (DEBUG) cat(file = stderr(), "scropius_trajectory_plot:NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scropius_trajectory_plot.RData", list = c(".schnappsEnv", ls(), ls(envir = globalenv())))
  }
  # cp =load(file="~/SCHNAPPsDebug/scropius_trajectory_plot.RData")
  
  vChanged = valuesChanged(parameters = c(
    "scorpDimX", "scorpDimY",
    "scorpMaxGenes", "scorpRepeat" , "scorpInFile"
  ))
  
  
  if (is.null(space) | vChanged) {
    if (vChanged) {
      cat(file = stderr(), "scropius Values changed\n")
    }
    require(ggplot2)
    p1 <- ggplot(projections, aes_string(dimX, dimY, colour = dimCol)) + geom_point()
    return(p1)
  }
  
  # space <- projections[, c(dimX, dimY)]
  require(SCORPIUS)
  # traj <- SCORPIUS::infer_trajectory(space)
  colnames(traj) <- c("Comp1", "Comp2", "time")
  draw_trajectory_plot(space, progression_group = projections[rownames(space), dimCol], path = as.matrix(traj[, 1:2]))
})

callModule(tableSelectionServer, "scorpiusTableMod", scorpiusModulesTable)
# selected clusters heatmap module

scorpiusHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusHeatmapPlotReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusHeatmapPlotReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scorpiusHeatmapPlotReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusHeatmapPlotReactive", id = "scorpiusHeatmapPlotReactive", duration = NULL)
  }
  
  # upI <- updateScorpiusInput() # needed to update input
  scEx = scEx()
  projections <- isolate(projections())
  traj <- scorpiusTrajectory()
  expr_sel <- scorpiusExpSel()
  modules <- scorpiusModules()
  
  dimCol <- input$dimScorpiusCol
  # doCalc <- input$scorpiusCalc
  pixelratio <- session$clientData$pixelratio
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  if (is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(traj)) {
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
  colnames(expr_sel$expr_sel) = rowData(scEx[colnames(expr_sel$expr_sel),])$symbol
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(x = expr_sel$expr_sel, time = traj$time, 
                                  progression_group = projections[rownames(expr_sel$expr_sel), dimCol], modules = modules,
                                  show_labels_row = TRUE, show_labels_col = FALSE,scale_features = TRUE,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {
    retVal
  })
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
    if (DEBUG) cat(file = stderr(), "downLoadTraj started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "downLoadTraj")
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "downLoadTraj")
      }
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("downLoadTraj", id = "downLoadTraj", duration = NULL)
    }
    
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


# Elpi Graph -----

# obeserver for elpi button ----
# observe: cellNameTable_rows_selected ----
observe(label = "ob_elpiCalcParams", {
  if (DEBUG) cat(file = stderr(), "observe elpiCalc\n")
  input$elpiCalc
  setRedGreenButtonCurrent(
    vars = list(
      c("elpinReps", input$elpinReps),
      c("elpiNumNodes", input$elpiNumNodes),
      c("elpiProbPoint", input$elpiProbPoint),
      c("ElpiMethod", input$ElpiMethod),
      c("dimElpi", input$dimElpi),
      c("elpiSeed", input$elpiSeed),
      c("dimElpi", input$dimElpi),
      c("dimElpiX", input$dimElpiX),
      c("dimElpiY", input$dimElpiY)
    )
  )
  
  updateButtonColor(buttonName = "elpiCalc", parameters = c(
    "elpinReps","elpiNumNodes","elpiProbPoint","ElpiMethod","dimElpi",
    "elpiSeed","dimElpi","dimElpiX","dimElpiY"
  ))
})

# elpiHeatmapPlotReactive ----
elpiHeatmapPlotReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "elpiHeatmapPlotReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiHeatmapPlotReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "elpiHeatmapPlotReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpiHeatmapPlotReactive", id = "elpiHeatmapPlotReactive", duration = NULL)
  }
  
  # upI <- updateScorpiusInput() # needed to update input
  scEx = scEx()
  clicked <- input$elpiCalc
  projections <- projections()
  psTime <- traj_getPseudotime()
  expr_sel <- traj_elpi_gimp()
  modules <- traj_elpi_modules()
  
  dimCol <- isolate(input$dimElpiCol)
  pixelratio <- session$clientData$pixelratio
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  
  if (is.null(projections) | is.null(modules) | is.null(expr_sel) | is.null(psTime)) {
    if (.schnappsEnv$DEBUG) cat(file = stderr(), paste("scorpiusHeatmapPlot:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/elpiHeatmapPlotReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/elpiHeatmapPlotReactive.RData")
  
  if (is.null(pixelratio)) pixelratio <- 1
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  colnames(expr_sel$expr_sel) = rowData(scEx[colnames(expr_sel$expr_sel),])$symbol
  
  outfile <- paste0(tempdir(), "/heatmapScorpius", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  outfile = NULL
  pst <- psTime$Pt[which(!is.na(psTime$Pt))]
  # modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
  retVal <- drawTrajectoryHeatmap(x =expr_sel$expr_sel,
                                  time = pst, progression_group = projections[rownames(expr_sel$expr_sel), dimCol], 
                                  modules = modules, show_labels_col = FALSE, show_labels_row = TRUE,
                                  filename = normalizePath(outfile, mustWork = FALSE)
  )
  
  exportTestValues(scorpiusHeatmapPlotReactive = {
    retVal
  })
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
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "elpi_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpi_plot", id = "elpi_plot", duration = NULL)
  }
  
  doCalc <- input$elpiCalc
  
  tree_data <- elpiTreeData()
  cep <- elpiGraphCompute()
  PointLabel <- elpiPointLabel()
  projections <- projections()
  dimX <- input$dimElpiX
  dimY <- input$dimElpiY
  dimCol <- input$dimElpiCol
  
  vChanged = valuesChanged(parameters = c(
    "elpinReps","elpiNumNodes","elpiProbPoint","ElpiMethod",
    "elpiSeed","dimElpi","dimElpiX","dimElpiY"
  ))
  
  
  if (is.null(projections)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpi_plot.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/elpi_plot.RData")
  if (is.null(cep) | vChanged) {
    if (vChanged) {
      cat(file = stderr(), "elpi Values changed\n")
    }
    require(ggplot2)
    p1 <- ggplot(projections, aes_string(dimX, dimY, colour = dimCol)) + geom_point()
    return(p1)
  }
  
  if (is.null(tree_data) | is.null(cep)) {
    return(NULL)
  }
  
  
  require(ggrepel)
  require(igraph)
  
  NodeLabs <- 1:nrow(cep[[length(cep)]]$NodePositions)
  NodeLabs[degree(ConstructGraph(cep[[length(cep)]])) != 1] <- NA
  
  p <- PlotPG(
    X = tree_data, TargetPG = cep[[length(cep)]],
    NodeLabels = NodeLabs,
    LabMult = 5, PointSize = NA, p.alpha = .5,
    GroupsLab = projections[rownames(tree_data),dimCol]
  )
  
  
  # p <- PlotPG(X = tree_data, TargetPG = cep[[length(cep)]], GroupsLab = PointLabel, p.alpha = 0.9)
  # p[[1]] <- p[[1]] + geom_label_repel(
  #   data = plyr::ddply(p[[1]]$data, ~Group, summarise, meanA = mean(PCA), meanB = mean(PCB)),
  #   aes(x = meanA, y = meanB, label = Group),
  #   vjust = 1
  # )
  p
})

# output$elpi_moduleHeatmap <- renderPlot({
#   draw_trajectory_heatmap(expr_sel, pst, group_name, modules)
# })

callModule(tableSelectionServer, "elpiTableMod", elpiModulesTable)

output$elpi_histo <- renderPlot({
  if (DEBUG) cat(file = stderr(), "elpi_histo started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpi_histo")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "elpi_histo")
    }
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
