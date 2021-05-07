
suppressMessages(require(Rsomoclu))
suppressMessages(require(kohonen))


# source("~/Rstudio/UTechSCB-SCHNAPPs/inst/app/serverFunctions.R")
# library(SCHNAPPs)
# library(shiny)

#' coE_somTrain
#' iData = expression matrix, rows = genes
#' cluster genes in SOM
#' returns genes associated together within a som-node
# coE_somTrain <- function(iData, nSom, geneName, ,clusterSOM = "dbCluster", clusterVal = "1") {
coE_somTrain <- function(iData, nSom, inputCells = "") {
  if (DEBUG) cat(file = stderr(), "coE_somTrain started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somTrain")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_somTrain")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somTrain", id = "coE_somTrain", duration = NULL)
  }
  
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somTrain.Rdata", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/coE_somTrain.Rdata")
  
  cols2use <- which (colnames(iData) %in% inputCells)
  res2 <- Rsomoclu.train(
    input_data = iData[, cols2use],
    nSomX = nSom, nSomY = nSom,
    nEpoch = 10,
    radius0 = 0,
    radiusN = 0,
    radiusCooling = "linear",
    mapType = "planar",
    gridType = "rectangular",
    scale0 = 1,
    scaleN = 0.01,
    scaleCooling = "linear"
  )
  colnames(res2$codebook) <- rownames(iData)[cols2use]
  rownames(res2$globalBmus) <- make.unique(as.character(rownames(iData)), sep = "___")
  return(res2)
} 
  
coE_somMap <- function(res2, iData) {
  if (DEBUG) cat(file = stderr(), "coE_somMap started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somMap")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_somMap")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somMap", id = "coE_somMap", duration = NULL)
  }
  if (is.null(res2) |is.null(iData)) {
    return(NULL)
  }
  sommap = Rsomoclu.kohonen(iData, res2)
  return(sommap)
}


# ## Show 'codebook'
# plot(sommap, type="codes", main = "Codes")
# ## Show 'component planes'
# plot(sommap, type = "property", property = sommap$codes[[1]][,1],
#      main = colnames(sommap$codes)[1])
# ## Show 'U-Matrix'
# plot(sommap, type="dist.neighbours")



coE_somGenes <- function(res2, geneName) {
  if (DEBUG) cat(file = stderr(), "coE_somGenes started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somGenes")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_somGenes")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somGenes", id = "coE_somGenes", duration = NULL)
  }
  if (is.null(res2) ) {
    return(NULL)
  }
  if (geneName %in% rownames(res2$globalBmus)) {
    return(NULL)
  }
  
  simGenes <- rownames(res2$globalBmus)[which(res2$globalBmus[, 1] == res2$globalBmus[geneName, 1] &
                                                res2$globalBmus[, 2] == res2$globalBmus[geneName, 2])]

  return(simGenes)
}

# SOM observers ----
.schnappsEnv$coE_SOMGrp <- "sampleNames"
observe(label = "ob13", {
  if (DEBUG) cat(file = stderr(), paste0("observe: coE_SOMGrp\n"))
  .schnappsEnv$coE_SOMGrp <- input$coE_SOMGrp
})
.schnappsEnv$coE_SOMSelection <- "1"
observe(label = "ob14", {
  if (DEBUG) cat(file = stderr(), paste0("observe: coE_SOMSelection\n"))
  .schnappsEnv$coE_SOMSelection <- input$coE_SOMSelection
})

# # coE_updateInputSOMt ====
# coE_updateInputSOMt <- reactive({
#   if (DEBUG) cat(file = stderr(), "coE_updateInputSOMt started.\n")
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "coE_updateInputSOMt")
#     if (!is.null(getDefaultReactiveDomain())) {
#       removeNotification(id = "coE_updateInputSOMt")
#     }
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("coE_updateInputSOMt", id = "coE_updateInputSOMt", duration = NULL)
#   }
#   
#   tsneData <- projections()
#   
#   # Can use character(0) to remove all choices
#   if (is.null(tsneData)) {
#     return(NULL)
#   }
#   
#   coln <- colnames(tsneData)
#   choices <- c()
#   for (cn in coln) {
#     if (length(levels(as.factor(tsneData[, cn]))) < 20) {
#       choices <- c(choices, cn)
#     }
#   }
#   if (length(choices) == 0) {
#     choices <- c("no valid columns")
#   }
#   updateSelectInput(
#     session,
#     "coE_clusterSOM",
#     choices = choices,
#     selected = .schnappsEnv$coE_SOMGrp
#   )
# })

# observeEvent(input$coE_clusterSOM,{
#   projections <- projections()
#   if (DEBUG) cat(file = stderr(), "observeEvent: input$coE_clusterSOM.\n")
#   # Can use character(0) to remove all choices
#   if (is.null(projections)) {
#     return(NULL)
#   }
#   if(!input$coE_clusterSOM %in% colnames(projections)) return(NULL)
#   choicesVal = levels(projections[, input$coE_clusterSOM])
#   updateSelectInput(
#     session,
#     "coE_clusterValSOM",
#     choices = choicesVal,
#     selected = .schnappsEnv$coE_SOMSelection
#   )
#   
# })



coE_somMapReact <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_somMapReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somMapReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_somMapReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somMapReact", id = "coE_somMapReact", duration = NULL)
  }
  
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmapWarning")
  }
  
  scEx_log <- scEx_log()
  input$updateSOMParameters
  
  selectedCells <- isolate(coE_SOM_dataInput())
  inputCells <- isolate(selectedCells$cellNames())
   res2 = isolate(coE_somTrainReact())
   
  if (is.null(scEx_log) | is.null(res2)) {
    return(NULL)
  }
   if (.schnappsEnv$DEBUGSAVE) {
     save(file = "~/SCHNAPPsDebug/coE_somMap.RData", list = c(ls()))
   }
   # cp = load(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData")
   
  iData <- as.matrix(assays(scEx_log)[[1]])
  cols2use <- which (colnames(iData) %in% inputCells)
  
  sommap <- coE_somMap(res2, iData[,cols2use]) 
  return(sommap)
})


coE_somTrainReact <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_somTrainReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somTrainReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_somTrainReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somTrainReact", id = "coE_somTrainReact", duration = NULL)
  }
  
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmapWarning")
  }
  
  scEx_log <- scEx_log()
  input$updateSOMParameters
  
  nSom <- isolate(input$coE_dimSOM)
  selectedCells <- isolate(coE_SOM_dataInput())
  inputCells <- isolate(selectedCells$cellNames())
  genesin <- isolate(input$coE_geneSOM)
  
  if (is.null(scEx_log)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somTrain.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData")
  
  iData <- as.matrix(assays(scEx_log)[[1]])
  cols2use <- which (colnames(iData) %in% inputCells)
  
  res2 = coE_somTrain(iData[,cols2use], nSom, inputCells)
  return(res2)
})



coE_somGenesReact <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_somGenesReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_somGenesReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_somGenesReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_somGenesReact", id = "coE_somGenesReact", duration = NULL)
  }
  
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmapWarning")
  }
  
  input$updateSOMParameters
  scEx_log <- scEx_log()
  
  genesin <- isolate(input$coE_geneSOM)
  res2 = isolate(coE_somTrainReact())
  
  if (is.null(scEx_log) | is.null(res2)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somGenes.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData")
  featureData <- rowData(scEx_log)
  geneName = geneName2Index(genesin, featureData)
  
  genes = coE_somGenes(res2, geneName)
  return(genes)
})

# coE_heatmapSOMReactive ----
#' coE_heatmapSOMReactive
#' calculates a self organizing map (SOM) on the expression data and identifies genes
#' that cluster together with a gene of interest
# TODO: check that we are using only raw counts and not normalized
coE_heatmapSOMReactive <- reactive({
  if (DEBUG) cat(file = stderr(), "coE_heatmapSOMReactive started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "coE_heatmapSOMReactive")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "coE_heatmapSOMReactive")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_heatmapSOMReactive", id = "coE_heatmapSOMReactive", duration = NULL)
  }
  
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmapWarning")
  }
  
  scEx_log <- scEx_log()
  projections <- projections()
  input$updateSOMParameters
  sampCol <- sampleCols$colPal
  ccols <- clusterCols$colPal
  # coE_updateInputSOMt()
  
  genesin <- isolate(input$coE_geneSOM)
  geneNames = isolate(coE_somGenesReact())
  selectedCells <- isolate(coE_SOM_dataInput())
  sampdesc <- isolate(selectedCells$selectionDescription())
  prj <- isolate(selectedCells$ProjectionUsed())
  prjVals <- isolate(selectedCells$ProjectionValsUsed())
  
  if (is.null(scEx_log) | is.null(geneNames) ) {
    return(
      list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "heatmap should be here"
      )
    )
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_heatmapSOMReactive.RData")
  
  scEx_matrix <- as.matrix(assays(scEx_log)[[1]])
  featureData <- rowData(scEx_log)
  # go from readable gene name to ENSG number
  genesin <- geneName2Index(genesin, featureData)
  

  # plot the genes found
  output$coE_somGenes <- renderText({
    paste(sampdesc, "\n",
          paste(featureData[geneNames, "symbol"], collapse = ", ", sep = ",")
    )
  })
  if (length(geneNames) < 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no additional gene found. reduce size of SOM",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }
  
  # create variables for heatmap module
  annotation <- data.frame(projections[, c(prj, "sampleNames")])
  rownames(annotation) <- rownames(projections)
  colnames(annotation) <- c(prj, "sampleNames")
  annCols <- list(
    "sampleNames" = sampCol,
    prj = ccols
  )
  
  retVal <- list(
    mat = scEx_matrix[geneNames, cellNs, drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # scale = "row",
    fontsize_row = 14,
    labels_row = featureData[geneNames, "symbol"],
    show_rownames = TRUE,
    annotation_col = annotation,
    show_colnames = FALSE,
    annotation_legend = TRUE,
    # breaks = seq(minBreak, maxBreak, by = stepBreak),
    # filename = 'test.png',
    # filename = normalizePath(outfile),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(
      n = 6, name =
        "RdBu"
    )))(6),
    annotation_colors = annCols
  )
  # system.time(do.call(TRONCO::pheatmap, retVal))
  
  setRedGreenButton(
    vars = list(
      c("coE_geneSOM", isolate(input$coE_geneSOM)),
      c("coE_dimSOM", nSOM),
      c("coE_SOM_dataInput-Mod_PPGrp", prjVals),
      c("coE_SOM_dataInput-Mod_clusterPP", prj)
    ),
    button = "updateSOMParameters"
  )
  
  exportTestValues(coE_heatmapSOMReactive = {
    retVal
  })
  return(retVal)
})
