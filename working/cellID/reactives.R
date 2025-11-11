# apply changes button 
# when clicked we actually do something
# observeEvent(input$updateCellID, )

library(CelliD)

cellIdReactive <- reactive({
  require(CelliD)
  if (DEBUG) cat(file = stderr(), "updateCellID started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "updateCellID")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "updateCellID")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("updateCellID", id = "updateCellID", duration = NULL)
  }
  
  scEx_log <- isolate(scEx_log())
  prjs <- isolate(sessionProjections$prjs)
  clicked <- input$updateCellID
  method <- isolate(input$cellID_Method)
  nmcs <- isolate(input$cellID_nmcs)
  
  if (is.null(scEx_log) ) {
    if (DEBUG) if (is.null(scEx_log)) cat(file = stderr(), "updateCellID scEx_log null.\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/updateCellID.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/updateCellID.RData")
  rownames(scEx_log) = make.names(make.unique(rowData(scEx_log)$symbol))
  switch(method,
         "MCA" = {
           scEx_log = RunMCA(scEx_log, nmcs = nmcs)
         }, 
         "DMAP" = {
           scEx_log = RunMCA(scEx_log, nmcs = nmcs)
           attr(reducedDim(scEx_log, "MCA"), "genesCoordinates") = attr(reducedDim(scEx_log, "MCA"), "genesCoordinates")[!duplicated(attr(reducedDim(scEx_log, "MCA"), "genesCoordinates"))]
           scEx_log = RunMCDMAP(X = scEx_log, dims = seq(nmcs), reduction.name = "MCA", k=10)
         }, 
         "TSNE" = {
           scEx_log = RunMCA(scEx_log, nmcs = nmcs)
           scEx_log = RunMCTSNE(scEx_log, dims =  seq(nmcs), reduction = "MCA", assay = "logcounts", reduction.name = "MCTSNE")
         }, 
         "UMAP" = {
           scEx_log = RunMCA(scEx_log, nmcs = nmcs)
           scEx_log = RunMCUMAP(scEx_log, dims =  seq(nmcs))
         }
  )
  # browser()
  rd = as.data.frame( reducedDim(scEx_log, "MCA")[colnames(scEx_log),])
  if (ncol(prjs) > 0) {
    prjs <- prjs[colnames(scEx_log), ]
    for (cn in colnames(rd)) {
      if (cn %in% colnames(prjs)) {
        prjs[, cn] <- rd[, cn]
      } else {
        prjs <- base::cbind(prjs, rd[, cn], deparse.level = 0)
        colnames(prjs)[ncol(prjs)] <- cn
      }
    }
    
  }else {
    prjs <- as.data.frame( rd)
  }
  # browser()
  sessionProjections$prjs <- prjs
  
  setRedGreenButton(
    vars = list(
      c("cellID_Method", method),
      c("cellID_nmcs", nmcs)
    ),
    button = "updateCellID"
  )
  
  .schnappsEnv$defaultValues[["cellID_Method"]] = isolate(input$cellID_Method)
  .schnappsEnv$defaultValues[["cellID_nmcs"]] = isolate(input$cellID_nmcs)
  
  
  # return the singleCellExperiment object with the redudced dimensions
  retVal <- scEx_log
  exportTestValues(cellIdReactive = {
    retVal
  })
  retVal
})

cellID_cellTableReact <- reactive({
  if (DEBUG) cat(file = stderr(), "cellID_cellTableReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "cellID_cellTableReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellID_cellTableReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("cellID_cellTableReact", id = "cellID_cellTableReact", duration = NULL)
  }
  mca = cellIdReactive()
  
  if (is.null(mca) ) {
    if (DEBUG) if (is.null(mca)) cat(file = stderr(), "cellID_cellTableReact mca null.\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/cellID_cellTableReact.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/cellID_cellTableReact.RData")
  
  retVal <- as.data.frame(reducedDim(mca, "MCA"))
  exportTestValues(cellID_cellTableReact = {
    retVal
  })
  retVal
})

cellID_geneTableReact <- reactive({
  if (DEBUG) cat(file = stderr(), "cellID_geneTableReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "cellID_geneTableReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellID_geneTableReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("cellID_geneTableReact", id = "cellID_geneTableReact", duration = NULL)
  }
  mca = cellIdReactive()
  
  if (is.null(mca) ) {
    if (DEBUG) if (is.null(mca)) cat(file = stderr(), "cellID_geneTableReact mca null.\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/cellID_geneTableReact.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/cellID_geneTableReact.RData")
  
  retVal <- as.data.frame(attr(reducedDim(mca, "MCA"), "genesCoordinates"))
  exportTestValues(cellID_geneTableReact = {
    retVal
  })
  retVal
  
})


cellID_stdevTableReact <- reactive({
  if (DEBUG) cat(file = stderr(), "cellID_stdevTableReact started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "cellID_stdevTableReact")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellID_stdevTableReact")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("cellID_stdevTableReact", id = "cellID_stdevTableReact", duration = NULL)
  }
  mca = cellIdReactive()
  
  if (is.null(mca) ) {
    if (DEBUG) if (is.null(mca)) cat(file = stderr(), "cellID_stdevTableReact mca null.\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/cellID_stdevTableReact.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/cellID_stdevTableReact.RData")
  
  retVal = data.frame(stdev = attr(reducedDim(mca, "MCA"), "stdev")) 
  exportTestValues(cellID_stdevTableReact = {
    retVal
  })
  retVal
})


# GetGroupGeneSet
cellIDGeneSetReact <- reactive({
  if (DEBUG) cat(file = stderr(), "cellIDGeneSet started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "cellIDGeneSet")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellIDGeneSet")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("cellIDGeneSet", id = "cellIDGeneSet", duration = NULL)
  }
  mca = cellIdReactive()
  gBy = input$cellID_gBy
  projections <- projections()
  
  if (is.null(mca) ) {
    if (DEBUG) if (is.null(mca)) cat(file = stderr(), "cellIDGeneSet mca null.\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/cellIDGeneSet.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/cellIDGeneSet.RData")
  # gBy = "sampleNames"
  colData(mca)[gBy] = projections[,gBy]
  # res = GetGroupGeneSet(mca, group.by = gBy, dims = seq(ncol(reducedDim(mca, "MCA"))))
  # GetCellGeneSet(mca, reduction = "MCA", cells = rownames(projections[projections$dbCluster=="1",]), n.features = 20, dims = c(2,3,4,5,6))
  res = GetGroupGeneRanking(mca, group.by = gBy, dims = seq(ncol(reducedDim(mca, "MCA"))))
  resFrame = data.frame(row.names = rownames(scEx_log))
  for (li in 1:length(res)) {
    resFrame[,names(res)[li]] = rep(NA, nrow(resFrame))
    resFrame[ names(res[[li]]),names(res[li])] = res[[li]]
  }
  # Just print text with the names
  retVal <- resFrame
  .schnappsEnv$defaultValues[["cellID_gBy"]] = isolate(input$cellID_gBy)
  
  exportTestValues(cellIDGeneSetReact = {
    retVal
  })
  retVal
  
})





