require(ElPiGraph.R)

# retVal <- drawTrajectoryHeatmap(x=expr_sel, time=traj$time, progression_group=projections[, dimCol], modules,
#                                 filename = normalizePath(outfile, mustWork = FALSE)
# )

# drawTrajectoryHeatmap ----
drawTrajectoryHeatmap <- function(x, time, progression_group = NULL, modules = NULL,
                                  show_labels_row = FALSE, show_labels_col = FALSE, scale_features = TRUE,
                                  ...) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(sQuote("x"), " must be a numeric matrix or data frame")
  }
  if (!is.vector(time) || !is.numeric(time)) {
    stop(sQuote("time"), " must be a numeric vector")
  }
  if (nrow(x) != length(time)) {
    stop(
      sQuote("time"), " must have one value for each row in ",
      sQuote("x")
    )
  }
  if ((!is.null(progression_group) && !is.vector(progression_group) &&
       !is.factor(progression_group)) || (!is.null(progression_group) &&
                                          length(progression_group) != nrow(x))) {
    stop(sQuote("progression_group"), " must be a vector or a factor of length nrow(x)")
  }
  if (is.null(rownames(x))) {
    rownames(x) <- paste("Row ", seq_len(nrow(x)))
  }
  col_ann <- data.frame(row.names = rownames(x), Time = time)
  x_part <- x[order(time), , drop = FALSE]
  if (scale_features) {
    x_part <- scale_quantile(x_part)
  }
  x_part <- t(x_part)
  gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }
  ann_col <- list(Time = RColorBrewer::brewer.pal(5, "RdGy"))
  if (!is.null(progression_group)) {
    if (!is.factor(progression_group)) {
      progression_group <- factor(progression_group)
    }
    col_ann$Progression <- progression_group
    num_progressions <- length(levels(progression_group))
    progression_cols <- if (num_progressions <= 9) {
      RColorBrewer::brewer.pal(num_progressions, "Set1")
    }
    else {
      gg_color_hue(num_progressions)
    }
    ann_col$Progression <- stats::setNames(
      progression_cols,
      levels(progression_group)
    )
  }
  labels_row <- if (!show_labels_row) {
    rep("", nrow(x_part))
  } else {
    NULL
  }
  labels_col <- if (!show_labels_col) {
    rep("", ncol(x_part))
  } else {
    NULL
  }
  if (!is.null(modules)) {
    x_part <- x_part[modules$feature, ]
    gaps_row <- which(modules$module[-1] != modules$module[-length(modules$module)])
    cluster_rows <- F
  }
  else {
    gaps_row <- NULL
    cluster_rows <- T
  }
  # return(list(data = x_part, cluster_cols = F, cluster_rows = cluster_rows,
  #        annotation_col = col_ann, annotation_colors = ann_col,
  #        gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col,
  #        ...))
  # pheatmap::pheatmap(x_part, cluster_cols = F, cluster_rows = cluster_rows,
  #                    annotation_col = col_ann, annotation_colors = ann_col,
  #                    gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col,
  #                    ...)
  list(
    mat = x_part, cluster_cols = F, cluster_rows = cluster_rows,
    annotation_col = col_ann, annotation_colors = ann_col,
    gaps_row = gaps_row, labels_row = labels_row, labels_col = labels_col
  )
}

# scorpiuseParameters ----
# herewith we control that scorpius is only run if the button is pressed.
scorpiuseParameters <- reactiveValues()

observeEvent(input$updatetScorpiusParameters,{
  # only react on this value
  cat(file = stderr(), paste("hit button\n"))
  scorpiuseParameters$dimX = isolate(input$dimScorpiusX)
  scorpiuseParameters$dimY = isolate(input$dimScorpiusY)
  scorpiuseParameters$scInput = isolate(scorpiusInput())
  scorpiuseParameters$scorpMaxGenes = isolate(input$scorpMaxGenes)
  scorpiuseParameters$scorpRepeat = isolate(input$scorpRepeat)
})


# scorpiusInput ----
# read input space from file
scorpiusInput <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusInput started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusInput")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusInput")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusInput", id = "scorpiusInput", duration = NULL)
  }
  
  inFile <- input$trajInputFile
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusInput.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # cp = load(file="~/SCHNAPPsDebug/scorpiusInput.RData")
  if (is.null(inFile)) {
    return(NULL)
  }
  if (!file.exists(inFile$datapath)) {
    return(NULL)
  }
  traj <- read.csv(file = inFile$datapath)
  if (colnames(traj) == c("path.Comp1", "path.Comp2", "time")) {
    return(traj)
  } else {
    warning("file not correct")
    return(NULL)
  }
})

# scorpiusSpace ----
# 2D space in which trajectory is calculated
# return projections or loaded space.
scorpiusSpace <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusSpace started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusSpace")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusSpace")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusSpace", id = "scorpiusSpace", duration = NULL)
  }
  
  
  projections <- isolate(projections())
  # doCalc <- input$scorpiusCalc
  dimX <- scorpiuseParameters$dimX
  dimY <- scorpiuseParameters$dimY
  scInput <- scorpiuseParameters$scInput
  
  if (!is.null(scInput)) {
    return(scInput[, c(1, 2)])
  }
  if ( is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusSpace:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusSpace.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusSpace.RData")
  
  space <- projections[, c(dimX, dimY)]
  
  return(space)
})


# scorpiusTrajectory ----
# calculate trajectory
scorpiusTrajectory <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusTrajectory started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusTrajectory")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusTrajectory")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusTrajectory", id = "scorpiusTrajectory", duration = NULL)
  }
  
  # only execute if button is pressed
  clicked = input$updatetScorpiusParameters
  cat(file = stderr(), paste("scorpiusTrajectory clicked:", clicked,"\n"))
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    return(NULL)
  }
  
  space <- isolate(scorpiusSpace())
  # doCalc <- input$scorpiusCalc
  scInput <- isolate(scorpiusInput())
  
  
  if (!is.null(scInput)) {
    return(scInput)
  }
  if (is.null(space)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusTrajectory:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusTrajectory.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusTrajectory.RData")
  require(SCORPIUS)
  traj <- SCORPIUS::infer_trajectory(space)
  traj$path = data.frame(traj$path)
  traj$path$idx <- 1:nrow(traj$path)
  traj$path <- traj$path[order(traj$time),] 
  rownames(traj$path) = names(sort(traj$time))
  traj$path$time = sort(traj$time)
  traj$path = traj$path[order(traj$path$idx),]
  traj$path = traj$path[, -3]
  # rownames(traj$path) == names(traj$time)
  
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
  
  return(traj$path)
})

# traj$path[1:10,]
# traj$path[order(traj$path$idx)[1:10],]

# scorpiusExpSel ----
scorpiusExpSel <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusExpSel started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusExpSel")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusExpSel")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusExpSel", id = "scorpiusExpSel", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification( id = "scorpiusExpSelWARNING")
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    return(NULL)
  }
  scEx_log <- isolate(scEx_log())
  traj <- isolate(scorpiusTrajectory())
  # doCalc <- input$scorpiusCalc
  scorpMaxGenes <- scorpiuseParameters$scorpMaxGenes
  scorpRepeat <- scorpiuseParameters$scorpRepeat
  
  if (is.null(traj) | is.null(scEx_log) | is.null(traj)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusExpSel:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusExpSel.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusExpSel.RData")
  cellsNotFound <- colnames(assays(scEx_log)[[1]])[!colnames(assays(scEx_log)[[1]]) %in% rownames(traj)]
  expression <- t(as.matrix(assays(scEx_log)[[1]][,rownames(traj)]))
  gimp <- gene_importances(expression[rownames(traj),], traj$time, num_permutations = scorpRepeat, num_threads = detectCores())
  maxRow <- min(scorpMaxGenes, nrow(gimp))
  gene_sel <- gimp[1:maxRow, ]
  expr_sel <- expression[, gene_sel$gene]
  # dfTmp = data.frame(matrix(0,nrow = length(cellsNotFound), ncol = ncol(expr_sel)))
  # rownames(dfTmp) = cellsNotFound
  # colnames(dfTmp) = colnames(expr_sel)
  # retVal = rbind(expr_sel,dfTmp)
  return(list(expr_sel = expr_sel, gene_sel = gene_sel))
})

# scorpiusModules ----
scorpiusModules <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusModules started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusModules")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusModules")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusModules", id = "scorpiusModules", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification( id = "scorpiusModulesWARNING")
  }
  
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    return(NULL)
  }

    # browser()
  scEx_log <- isolate(scEx_log())
  # projections = projections()
  # space <- scorpiusSpace()
  traj <- isolate(scorpiusTrajectory())
  expr_sel <- isolate(scorpiusExpSel())
  
  # scorpiusModules = scorpiusModules()
  # upI <- updateScorpiusInput() # needed to update input
  dimX <- scorpiuseParameters$dimX
  dimY <- scorpiuseParameters$dimY
  # dimCol = input$dimScorpiusCol
  # doCalc <- input$scorpiusCalc
  
  if (is.null(traj) | is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusModules:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusModules.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusModules.RData")
  
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(scEx_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  
  modules <- extract_modules(scale_quantile(expr_sel$expr_sel), traj$time, verbose = T)
  modules <- as.data.frame(modules)
  fd <- rowData(scEx_log)
  modules$symbol <- fd[modules$feature, "symbol"]
  rownames(modules) <- make.unique(as.character(modules$symbol, sep = "___"))
  return(modules)
})

# scorpiusModulesTable ----
scorpiusModulesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "scorpiusModulesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "scorpiusModulesTable")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scorpiusModulesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scorpiusModulesTable", id = "scorpiusModulesTable", duration = NULL)
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetScorpiusParameters == 0) {
    return(NULL)
  }
  
  scEx_log <- isolate(scEx_log())
  # projections = projections()
  # space <- scorpiusSpace()
  traj <- isolate(scorpiusTrajectory())
  expr_sel <- isolate(scorpiusExpSel())
  modules <- isolate(scorpiusModules())
  
  # scorpiusModules = scorpiusModules()
  # upI <- updateScorpiusInput() # needed to update input
  dimX <- scorpiuseParameters$dimX
  dimY <- scorpiuseParameters$dimY
  # dimCol = input$dimScorpiusCol
  # doCalc <- input$scorpiusCalc
  
  if (is.null(traj) | is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusModules:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusModulesTable.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusModulesTable.RData")
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(scEx_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  
  gene_selDF <- as.data.frame(expr_sel$gene_sel)
  rownames(gene_selDF) = gene_selDF[,1]
  gene_selDF = gene_selDF[,-1]
  return(cbind(modules,gene_selDF[modules$feature,]))
})



# --------------------------
# Elpi Graph
# --------------------------
# TODO change to projections input


# elpiTreeData ----
elpiTreeData <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiTreeData\n")
  }
  
  scEx <- scEx()
  projections <- projections()
  dimElpi <- isolate(input$dimElpi)
  dim1 <- input$dimElpiX
  dim2 <- input$dimElpiY
  if (is.null(scEx)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiTreeData.RData", list = c(base::ls()))
  }
  # cp =load(file = "~/SCHNAPPsDebug/elpiTreeData.RData")
  if (dimElpi == "elpiPCA") {
    return(t(as.matrix(assays(scEx)[[1]])))
  }
  if (dimElpi == "components") {
    return(as.matrix(projections[,c(dim1,dim2)]))
  }
  cat(file = stderr(), "elpiTreeData should not happen\n")
  t(assays(scEx)[[1]])
})


# traj_endpoints ----
traj_endpoints <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_endpoints")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_endpoints")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_endpoints", id = "traj_endpoints", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_endpoints started.\n")
  scEx_log <- scEx_log()
  projections <- projections()
  TreeEPG <- elpiGraphCompute()
  elpimode <- input$ElpiMethod
  seed <- input$elpiSeed
  
  if (is.null(projections) || is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle") {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_endpoints.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/traj_endpoints.RData")
  set.seed(seed = seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[length(TreeEPG)]])
  Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end', Circular = T,Nodes = 30,KeepEnds = T)
  # get all end-points:
  endPoints = unique(c(sapply(Tree_e2e, function(x) x[1]), sapply(Tree_e2e, function(x) x[length(x)])))
  return(endPoints)
})


# updateElpiInput ----
# update the start/end positions with the possible end point values
# traj_getPseudotime ----
traj_getPseudotime <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_getPseudotime")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_getPseudotime")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_getPseudotime", id = "traj_getPseudotime", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
  scEx_log <- scEx_log()
  # scEx_log_sha <- scEx_log_sha()
  TreeEPG <- elpiGraphCompute()
  TreeEPG_sha <- TreeEPG_sha()
  
  elpimode <- input$ElpiMethod
  tree_data <- elpiTreeData()
  tragetPath <- traj_tragetPath()
  targetPathSha <- traj_targetPathSha()
  seed <- input$elpiSeed
  
  if (is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle" || length(tragetPath) == 0) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_getPseudotime.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/traj_getPseudotime.RData")
  
  # cacheResult <- checkShaCache(moduleName = "traj_getPseudotime", 
  #                              moduleParameters = list(scEx_log_sha, TreeEPG_sha,
  #                                                      tree_data, targetPathSha,
  #                                                      elpimode, seed))
  # if (cacheResult$status == "finished") {
  #   return(cacheResult$retVal)
  # }
  # if (cacheResult$status == "running") {
  #   showNotification(cacheResult$message, id = "traj_elpi_modules", duration = 10)
  #   return(NULL)
  # }
  # if (cacheResult$status == "error") {
  #   showNotification(cacheResult$message, id = "traj_elpi_modules", duration = NULL, type="error")
  #   return(NULL)
  # }
  
  set.seed(seed)
  # computing a Partition structure
  PartStruct <- ElPiGraph.R::PartitionData(X = tree_data, NodePositions = TreeEPG[[length(TreeEPG)]]$NodePositions)
  
  # projection structure
  ProjStruct <- ElPiGraph.R::project_point_onto_graph(X = tree_data,
                                                      NodePositions = TreeEPG[[length(TreeEPG)]]$NodePositions,
                                                      Edges = TreeEPG[[length(TreeEPG)]]$Edges$Edges,
                                                      Partition = PartStruct$Partition)
  psTime = ElPiGraph.R::getPseudotime(ProjStruct = ProjStruct, NodeSeq = names(tragetPath[[1]]))
  
  # writeShaCache(moduleName = "traj_getPseudotime", 
  #               moduleParameters = list(scEx_log_sha, TreeEPG,
  #                                       tree_data, targetPathSha,
  #                                       elpimode, seed),
  #               retVal = psTime,
  #               status = "finished",
  #               message = "")
  return(psTime)
  
})

traj_targetPathSha <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_targetPathSha")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_targetPathSha")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_targetPathSha", id = "traj_targetPathSha", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
  # browser()
  targetPath <- traj_tragetPath()
  if(is.null(targetPath) | length(targetPath)==0)
    return("")
  sha1(as.character(targetPath[[1]]))
})

## traj_elpi_modules -----
traj_elpi_modules <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_elpi_modules")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_elpi_modules")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_elpi_modules", id = "traj_elpi_modules", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_elpi_modules started.\n")
  scEx_log <- scEx_log()
  # scEx_log_sha<- scEx_log_sha()
  TreeEPG <- elpiGraphCompute()
  TreeEPG_sha <- TreeEPG_sha()
  
  elpimode <- input$ElpiMethod
  gene_sel <- traj_elpi_gimp()
  seed <- input$elpiSeed
  
  # add a button to check result.
  # This button will invalidate this reactive and restart checking of cache.
  # input$elpi_modules_check
  
  if (is.null(scEx_log) | is.null(gene_sel) | elpimode=="computeElasticPrincipalCircle" ) {
    return(NULL)
  }
  if (nrow(gene_sel$gene_sel)<1) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_elpi_modules.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/traj_elpi_modules.RData")
  
  # cacheResult <- checkShaCache(moduleName = "traj_elpi_modules", 
  #                             moduleParameters = list(scEx_log_sha, TreeEPG_sha,
  #                                                     elpimode, gene_sel, seed))
  # if (cacheResult$status == "finished") {
  #   return(chacheResult$retVal)
  # }
  # if (cacheResult$status == "running") {
  #   showNotification(cacheResult$message, id = "traj_elpi_modules", duration = 10)
  #   return(NULL)
  # }
  # if (cacheResult$status == "error") {
  #   showNotification(cacheResult$message, id = "traj_elpi_modules", duration = NULL, type="error")
  #   return(NULL)
  # }
  gene_sel = gene_sel$gene_sel
  set.seed(seed)
  expr_sel <- t(as.matrix(assays(scEx_log)[[1]][gene_sel$gene,]))
  
  ## Group the genes into modules and visualise the modules in a heatmap
  # group_name should be dbCluster or other selectable option
  modules <- SCORPIUS::extract_modules(SCORPIUS::scale_quantile(expr_sel))
  modules <- as.data.frame(modules)
  fd <- rowData(scEx_log)
  modules$symbol <- fd[modules$feature, "symbol"]
  rownames(modules) <- make.unique(as.character(modules$symbol, sep = "___"))
  
  # writeShaCache(moduleName = "traj_elpi_modules", 
  #               moduleParameters = list(scEx_log_sha, TreeEPG_sha,
  #                                       elpimode, gene_sel, seed),
  #               retVal = modules,
  #               status = "finished",
  #               message = "")
  return(modules)
})

TreeEPG_sha <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "TreeEPG_sha")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "TreeEPG_sha")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("TreeEPG_sha", id = "TreeEPG_sha", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "TreeEPG_sha started.\n")
  # browser()
  TreeEPG <- elpiGraphCompute()
  if (is.null(TreeEPG)){
    return(NULL)
  }
  sha1(as.character(TreeEPG))
})

# traj_elpi_gimp -----
traj_elpi_gimp <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_elpi_gimp")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_elpi_gimp")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_elpi_gimp", id = "traj_elpi_gimp", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_elpi_gimp started.\n")
  scEx_log <- scEx_log()
  projections <- projections()
  TreeEPG <- elpiGraphCompute()
  elpimode <- input$ElpiMethod
  psTime = traj_getPseudotime()
  num_permutations <- input$elpi_num_permutations
  ntree <- input$elpi_ntree
  ntree_perm <- input$elpi_ntree_perm
  nGenes <- input$elpi_nGenes
  seed <- input$elpiSeed
  
  if (is.null(scEx_log) || is.null(psTime) || elpimode=="computeElasticPrincipalCircle") {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_elpi_gimp.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/traj_elpi_gimp.RData")
  
  set.seed(seed)
  require(SCORPIUS)
  logCounts <- as.matrix(assays(scEx_log)[[1]][,which(!is.na(psTime$Pt))])
  pst = psTime$Pt[which(!is.na(psTime$Pt))]
  geneImport <- SCORPIUS::gene_importances(t(logCounts), pst, num_permutations = num_permutations, ntree = ntree,
                                           ntree_perm = ntree_perm, mtry = ncol(logCounts) * 0.01, num_threads = detectCores()-1)
  gene_sel <- geneImport[1:nGenes,]
  expr_sel <- t(logCounts)[, gene_sel$gene]
  
  return(list(expr_sel = expr_sel, gene_sel = gene_sel))
})
# 
# # observeProj ----
# # update projections
# observe({
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "observeProj")
#     if (!is.null(getDefaultReactiveDomain()))
#       removeNotification(id = "observeProj")
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("observeProj", id = "observeProj", duration = NULL)
#   }
#   if (DEBUG) cat(file = stderr(), "observeProj started.\n")
#   
#   startNode <- input$elpiStartNode
#   endNode <- input$elpiEndNode
#   elpimode <- input$ElpiMethod
#   psTime = traj_getPseudotime()
#   scEx_log <- scEx_log()
#   isolate({
#     prjs <- sessionProjections$prjs
#   })
#   
#   if (is.null(scEx_log) || is.null(psTime) || elpimode=="computeElasticPrincipalCircle") {
#     return(NULL)
#   }
#   if (.schnappsEnv$DEBUGSAVE) {
#     save(file = "~/SCHNAPPsDebug/observeProj.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/SCHNAPPsDebug/observeProj.RData")
#   cn = paste0("traj_", startNode, "_", endNode)
#   if (cn %in% colnames(prjs)) {
#     return(NULL)
#   }
#   # browser()
#   if (ncol(prjs) > 0) {
#     # make sure we are working with the correct cells. This might change when cells were removed.
#     prjs = prjs[colnames(scEx_log),,drop=FALSE]
#     # didn't find a way to easily overwrite columns
#     
#     if (cn %in% colnames(prjs)) {
#       prjs[, cn] <- psTime$Pt
#     } else {
#       prjs <- base::cbind(prjs, psTime$Pt, deparse.level = 0)
#       colnames(prjs)[ncol(prjs)] <- cn
#     }
#     sessionProjections$prjs <- prjs
#   } else {
# 
#         prjs <- data.frame(cn = psTime$Pt )
#     rownames(prjs) = colnames(scEx_log)
#     colnames(prjs)[ncol(prjs)] = cn
#     sessionProjections$prjs = prjs
#   }
# })

#traj_tragetPath ----
traj_tragetPath <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "traj_tragetPath")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "traj_tragetPath")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("traj_tragetPath", id = "traj_tragetPath", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "traj_tragetPath started.\n")
  scEx_log <- scEx_log()
  projections <- projections()
  TreeEPG <- elpiGraphCompute()
  elpimode <- input$ElpiMethod
  startNode <- input$elpiStartNode
  endNode <- input$elpiEndNode
  seed <- input$elpiSeed
  
  if (is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle") {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_tragetPath.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/traj_tragetPath.RData")
  set.seed(seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[length(TreeEPG)]])
  Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end')
  
  retVal <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[1] %in% c(startNode,endNode)) & any(x[length(x)] %in% c(startNode,endNode)) })]
  
})

# elpiGraphCompute ----
elpiGraphCompute <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiGraphCompute\n")
  }
  
  doCalc <- input$elpiCalc
  tree_data <- isolate(elpiTreeData())

  nReps <- isolate(input$elpinReps) # 1-50
  NumNodes <- isolate(input$elpiNumNodes) # 10 - 100
  ProbPoint <- isolate(input$elpiProbPoint) # 0.1-1.0
  method <- isolate(input$ElpiMethod)
  dimUse <-  isolate(input$dimElpi)
  seed <- isolate(input$elpiSeed)
  
  require(parallel)
  
  if (is.null(tree_data)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiCalc.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load("~/SCHNAPPsDebug/elpiCalc.RData")
  
  #   cacheResult <- checkShaCache(moduleName = "elpiGraphCompute", 
  #                              moduleParameters = list(nReps, NumNodes, ProbPoint,
  #                                                      method, doCalc, dimUse, seed, tree_data))
  # if (cacheResult$status == "finished") {
  #   return(cacheResult$retVal)
  # }
  # if (cacheResult$status == "running") {
  #   showNotification(cacheResult$message, id = "traj_elpi_modules", duration = 10)
  #   return(NULL)
  # }
  # if (cacheResult$status == "error") {
  #   showNotification(cacheResult$message, id = "traj_elpi_modules", duration = NULL, type="error")
  #   return(NULL)
  # }
  
  set.seed(seed)
  if (dimUse == "elpiPCA") {
    elpiDoPCA = TRUE
    elipCenter = TRUE
  } else {
    elpiDoPCA = FALSE
    elipCenter = FALSE
    
  }
  cep <- do.call(method, list(
    X = tree_data,
    NumNodes = NumNodes,
    drawAccuracyComplexity = F,
    nReps = nReps, # bootstrapping
    ProbPoint = ProbPoint, # bootstrapping
    drawPCAView = F,
    drawEnergy = F,
    Do_PCA = elpiDoPCA,
    CenterData = elipCenter, 
    n.cores = detectCores() - 1
  ))
  # writeShaCache(moduleName = "elpiGraphCompute", 
  #               moduleParameters = list(nReps, NumNodes, ProbPoint,
  #                                       method, doCalc, dimUse, seed, tree_data),
  #               retVal = cep,
  #               status = "finished",
  #               message = "")
  setRedGreenButton(
    vars = list(
      c("elpinReps", isolate(input$elpinReps)),
      c("elpiNumNodes", isolate(input$elpiNumNodes)),
      c("elpiProbPoint", isolate(input$elpiProbPoint)),
      c("ElpiMethod", isolate(input$ElpiMethod)),
      c("dimElpi", isolate(input$dimElpi)),
      c("elpiSeed", isolate(input$elpiSeed)),
      c("dimElpi", isolate(input$dimElpi)),
      c("dimElpiX", isolate(input$dimElpiX)),
      c("dimElpiY", isolate(input$dimElpiY))
    ),
    button = "elpiCalc"
  )
  
  return(cep)
})

# elpiGraphConstruct ----
elpiGraphConstruct <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiGraphConstruct\n")
  }
  cep <- elpiGraphCompute()
  tree_data <- elpiTreeData()
  seed <- input$elpiSeed
  
  if (is.null(cep)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiConstruct.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/elpiConstruct.RData")
  
  set.seed(seed = seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(PrintGraph = cep[[length(cep)]])
  Tree_Brches <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = "branches")
  PartStruct <- ElPiGraph.R::PartitionData(X = tree_data, NodePositions = cep[[length(cep)]]$NodePositions)
  
  list(
    Tree_Graph = Tree_Graph,
    Tree_e2e = GetSubGraph(Net = Tree_Graph, Structure = "end2end"),
    Tree_Brches = Tree_Brches,
    Tree_BrBrPt = GetSubGraph(Net = Tree_Graph, Structure = "branches&bpoints"),
    PartStruct = PartStruct,
    PtInBr = lapply(Tree_Brches, function(x) {
      which(PartStruct$Partition %in% x)
    })
  )
})

# elpiPointLabel ----
elpiPointLabel <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "elpiPointLabel\n")
  }
  elpiGraphConstruct <- elpiGraphConstruct()
  if (is.null(elpiGraphConstruct)) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiPointLabel.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load(file = "~/SCHNAPPsDebug/elpiPointLabel.RData")
  PartStruct <- elpiGraphConstruct$PartStruct
  Tree_BrBrPt <- elpiGraphConstruct$Tree_BrBrPt
  
  
  PointLabel <- rep("", length(PartStruct$Partition))
  
  for (i in 1:length(Tree_BrBrPt)) {
    PointLabel[PartStruct$Partition %in% Tree_BrBrPt[[i]]] <- names(Tree_BrBrPt)[i]
  }
  return(PointLabel)
})

# elpiModulesTable ----
elpiModulesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "elpiModulesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiModulesTable")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpiModulesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("elpiModulesTable", id = "elpiModulesTable", duration = NULL)
  }
  scEx_log <- scEx_log()
  # projections = projections()
  # space <- scorpiusSpace()
  traj <- traj_getPseudotime()
  expr_sel <- traj_elpi_gimp()
  modules <- traj_elpi_modules()
  
  # scorpiusModules = scorpiusModules()
  # upI <- updateScorpiusInput() # needed to update input
  dimX <- input$dimElpiX
  dimY <- input$dimElpiY
  # dimCol = input$dimScorpiusCol
  doCalc <- input$elpiCalc
  
  if (!doCalc | is.null(scEx_log) | is.null(expr_sel)) {
    if (DEBUG) cat(file = stderr(), paste("elpiModulesTable:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/elpiModulesTable.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/elpiModulesTable.RData")
  # space = projections[,c(dimX, dimY)]
  # traj <- infer_trajectory(space)
  # expression = as.matrix(exprs(scEx_log))
  # gimp <- gene_importances(t(expression), traj$time, num_permutations = 0, num_threads = 8)
  # gene_sel <- gimp[1:50,]
  # expr_sel <- t(expression)[,gene_sel$gene]
  
  # browser()
  gene_selDF <- as.data.frame(expr_sel$gene_sel)
  rownames(gene_selDF) = gene_selDF[,1]
  gene_selDF = gene_selDF[,-1]
  return(cbind(modules,gene_selDF[modules$feature,]))
})


# here we define reactive values/variables
# e.g.
# DummyReactive = reactive({
#   # some debugging messages
#   if(DEBUG)cat(file=stderr(), "pca\n")
#   # call dependancies (reactives)
#   gbm_log = gbm_log()
#   # check if they are available
#   if(is.null(gbm_log)){
#     if(DEBUG)cat(file=stderr(), "pca:NULL\n")
#     return(NULL)
#   }
#   if(!is.null(getDefaultReactiveDomain())){
#     showNotification("loading", id="DummyFunc", duration = NULL)
#   }
#   if(DEBUGSAVE)
#     save(file = "~/SCHNAPPsDebug/DummyReactive.RData", list = c(ls(),ls(envir = globalenv())))
#   # load(file='~/SCHNAPPsDebug/DummyReactive.RData')
#
#   # actual calculation
#   retVal = DummyFunc(gbm_log)
#
#   if(retVal == 0 & !is.null(getDefaultReactiveDomain()) ){
#     showNotification("Dummy is 0", type="warning", duration=NULL) # has to be removed by use, no removeNotification is following.
#     return(NULL)
#   }
#   if(DEBUG)cat(file=stderr(), "inputData: done\n")
#   if(!is.null(getDefaultReactiveDomain())){
#     removeNotification( id="DummyFunc")
#   }
#   return(retVal)
# })
#
# # declare heavy calculations
# myHeavyCalculations = list(c("running DummyReactive", "DummyReactive"))


# this will never be executed as it won't be called...

# imageDummyPrecompute <- reactive({
#   if(DEBUG)cat(file=stderr(), "imageDummyPrecompute\n")
#
#   gbm = gbm()
#   # check if they are available
#   if(is.null(gbm)){
#     if(DEBUG)cat(file=stderr(), "imageDummyPrecompute:NULL\n")
#     return(NULL)
#   }
#   width  <- session$clientData$output_plot_width
#   height <- session$clientData$output_plot_height
#
#
#   if(!is.null(getDefaultReactiveDomain())){
#     showNotification("loading", id="DummyFunc", duration = NULL)
#   }
#   # For high-res displays, this will be greater than 1
#   pixelratio <- session$clientData$pixelratio
#   if(is.null(pixelratio)) pixelratio = 1
#   width  <- session$clientData$output_plot_width
#   height <- session$clientData$output_plot_height
#   if(is.null(width)){width=96*7} # 7x7 inch output
#   if(is.null(height)){height=96*7}
#
#   myPNGwidth <- width/96
#   myPNGheight <- height/96
#
#   outfile <- paste0(tempdir(),'/dummy.png')
#   if(DEBUG)cat(file=stderr(), paste("output file: ", outfile, "\n"))
#   if(DEBUG)cat(file=stderr(), paste("output file normalized: ", normalizePath(outfile), "\n"))
#   m = data.frame("V1"=colSums(as.matrix(exprs(gbm))))
#   p=ggplot(m, aes(V1)) + geom_bar()
#   ggsave(file=normalizePath(outfile), plot=p, width=myPNGwidth, height=myPNGheight, units="in")
#
#   if(DEBUG)cat(file=stderr(), "imageDummyPrecompute:done\n")
#
#   return(list(src = normalizePath(outfile),
#               contentType = 'image/png',
#               width = width,
#               height = height,
#               alt = "Dummy should be here"))
#
#
# })
