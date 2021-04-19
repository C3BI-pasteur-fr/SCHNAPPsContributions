require(ElPiGraph.R)
require(Tempora)
require(S4Vectors)
require(SingleCellExperiment)

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
    x_part <- x_part[modules$symbol, ]
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
    save(file = "~/SCHNAPPsDebug/scorpiusInput.RData", list = c(ls()))
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
  
  
  # doCalc <- input$scorpiusCalc
  dimX <- scorpiuseParameters$dimX
  dimY <- scorpiuseParameters$dimY
  scInput <- scorpiuseParameters$scInput
  projections <- projections()
  
  if (!is.null(scInput)) {
    return(scInput[, c(1, 2)])
  }
  if ( is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusSpace:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusSpace.RData", list = c(ls()))
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
  
  space <- scorpiusSpace()
  # doCalc <- input$scorpiusCalc
  scInput <- scorpiusInput()
  
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/scorpiusTrajectory.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/scorpiusTrajectory.RData")
  if (!is.null(scInput)) {
    return(scInput)
  }
  if (is.null(space)) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusTrajectory:NULL\n"))
    return(NULL)
  }
  if (nrow(space)<10 | ncol(space)<2) {
    if (DEBUG) cat(file = stderr(), paste("scorpiusTrajectory:NULL; need more samples/columns\n"))
    return(NULL)
  }
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
    save(file = "~/SCHNAPPsDebug/scorpiusExpSel.RData", list = c(ls()))
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
    save(file = "~/SCHNAPPsDebug/scorpiusModules.RData", list = c(ls()))
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
    save(file = "~/SCHNAPPsDebug/scorpiusModulesTable.RData", list = c(ls()))
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
  
  clicked <- input$elpiCalc
  scEx <- scEx()
  projections <- isolate(projections())
  dimElpi <- isolate(input$dimElpi)
  dim1 <- input$dimElpiX
  dim2 <- input$dimElpiY
  if (is.null(scEx)) {
    cat(file = stderr(), "--- elpiTreeData: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiTreeData.RData", list = c(base::ls()))
  }
  # cp =load(file = "~/SCHNAPPsDebug/elpiTreeData.RData")
  if (dimElpi == "elpiPCA") {
    return(t(as.matrix(assays(scEx)[[1]])))
  }
  if (!all(c(dim1, dim2) %in% colnames(projections))){
    cat(file = stderr(), paste("--- dims not in projctions: NULL:", c(dim1, dim2), "\n"))
    return(NULL)
    
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
  
  clicked <- input$elpiCalc
  TreeEPG <- elpiGraphCompute()
  scEx_log <- scEx_log()
  
  projections <- isolate(projections())
  elpimode <- isolate(input$ElpiMethod)
  seed <- isolate(input$elpiSeed)
  
  if (is.null(projections) || is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle") {
    cat(file = stderr(), "--- traj_endpoints: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_endpoints.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/traj_endpoints.RData")
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
  clicked <- input$elpiCalc
  scEx_log <- scEx_log()
  # scEx_log_sha <- scEx_log_sha()
  TreeEPG <- elpiGraphCompute()
  # TreeEPG_sha <- TreeEPG_sha()
  
  tree_data <- elpiTreeData()
  tragetPath <- traj_tragetPath()
  seed <- isolate(input$elpiSeed)
  elpimode <- isolate(input$ElpiMethod)
  
  if (is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle" || length(tragetPath) == 0) {
    cat(file = stderr(), "--- traj_getPseudotime: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_getPseudotime.RData", list = c(ls()))
  }
  # cp =load(file="~/SCHNAPPsDebug/good/traj_getPseudotime.RData")
  
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

# traj_targetPathSha <- reactive({
#   start.time <- base::Sys.time()
#   on.exit({
#     printTimeEnd(start.time, "traj_targetPathSha")
#     if (!is.null(getDefaultReactiveDomain()))
#       removeNotification(id = "traj_targetPathSha")
#   })
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("traj_targetPathSha", id = "traj_targetPathSha", duration = NULL)
#   }
#   if (DEBUG) cat(file = stderr(), "traj_getPseudotime started.\n")
#   # browser()
#   targetPath <- traj_tragetPath()
#   if(is.null(targetPath) | length(targetPath)==0)
#     return("")
#   sha1(as.character(targetPath[[1]]))
# })

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
  
  clicked <- input$elpiCalc
  scEx_log <- scEx_log()
  TreeEPG <- elpiGraphCompute()
  gene_sel <- traj_elpi_gimp()
  
  elpimode <- isolate(input$ElpiMethod)
  seed <- isolate(input$elpiSeed)
  
  # add a button to check result.
  # This button will invalidate this reactive and restart checking of cache.
  # input$elpi_modules_check
  
  if (is.null(scEx_log) | is.null(gene_sel) | elpimode=="computeElasticPrincipalCircle" ) {
    cat(file = stderr(), "--- traj_elpi_modules: NULL\n")
    return(NULL)
  }
  if (nrow(gene_sel$gene_sel)<1) {
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_elpi_modules.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/traj_elpi_modules.RData")
  
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
  
  return(modules)
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
  psTime = traj_getPseudotime()
  clicked <- input$elpiCalc
  
  elpimode <- isolate(input$ElpiMethod)
  num_permutations <- input$elpi_num_permutations
  ntree <- input$elpi_ntree
  ntree_perm <- input$elpi_ntree_perm
  nGenes <- input$elpi_nGenes
  seed <- isolate(input$elpiSeed)
  
  if (is.null(scEx_log) || is.null(psTime) || elpimode=="computeElasticPrincipalCircle") {
    cat(file = stderr(), "--- traj_elpi_gimp: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_elpi_gimp.RData", list = c(ls()))
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
#     save(file = "~/SCHNAPPsDebug/observeProj.RData", list = c(ls()))
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
  clicked <- input$elpiCalc
  startNode <- input$elpiStartNode
  endNode <- input$elpiEndNode
  
  elpimode <- isolate(input$ElpiMethod)
  seed <- isolate(input$elpiSeed)
  
  if (is.null(scEx_log) || is.null(TreeEPG) || elpimode=="computeElasticPrincipalCircle") {
    cat(file = stderr(), "--- traj_tragetPath: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/traj_tragetPath.RData", list = c(ls()))
  }
  # cp = load(file="~/SCHNAPPsDebug/traj_tragetPath.RData")
  set.seed(seed)
  Tree_Graph <- ElPiGraph.R::ConstructGraph(TreeEPG[[length(TreeEPG)]])
  Tree_e2e <- ElPiGraph.R::GetSubGraph(Net = Tree_Graph, Structure = 'end2end')
  
  retVal <- Tree_e2e[sapply(Tree_e2e, function(x){any(x[1] %in% c(startNode,endNode)) & any(x[length(x)] %in% c(startNode,endNode)) })]
  
})

# elpiGraphCompute ----
elpiGraphCompute <- reactive({
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "elpiGraphCompute")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "elpiGraphCompute")
  })
  if (DEBUG) {
    cat(file = stderr(), "elpiGraphCompute\n")
  }
  
  clicked <- input$elpiCalc
  tree_data <- isolate(elpiTreeData())
  
  nReps <- isolate(input$elpinReps) # 1-50
  NumNodes <- isolate(input$elpiNumNodes) # 10 - 100
  ProbPoint <- isolate(input$elpiProbPoint) # 0.1-1.0
  method <- isolate(input$ElpiMethod)
  dimUse <-  isolate(input$dimElpi)
  seed <- isolate(input$elpiSeed)
  
  require(parallel)
  
  if (is.null(tree_data)) {
    cat(file = stderr(), "--- elpiGraphCompute: NULL\n")
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    base::save(file = "~/SCHNAPPsDebug/elpiCalc.RData", list = c(base::ls(), base::ls(envir = globalenv())))
  }
  # load("~/SCHNAPPsDebug/elpiCalc.RData")
  
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
  clicked <- input$elpiCalc
  cep <- elpiGraphCompute()
  tree_data <- elpiTreeData()
  
  seed <- isolate(input$elpiSeed)
  
  if (is.null(cep)) {
    cat(file = stderr(), "--- elpiGraphConstruct: NULL\n")
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
  
  clicked <- input$elpiCalc
  elpiGraphConstruct <- elpiGraphConstruct()
  if (is.null(elpiGraphConstruct)) {
    cat(file = stderr(), "--- elpiPointLabel: NULL\n")
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
  clicked <- input$elpiCalc
  
  scEx_log <- isolate(scEx_log())
  traj <- traj_getPseudotime()
  expr_sel <- traj_elpi_gimp()
  modules <- traj_elpi_modules()
  
  dimX <- isolate(input$dimElpiX)
  dimY <- isolate(input$dimElpiY)
  
  if (is.null(scEx_log) | is.null(expr_sel)) {
    if (DEBUG) cat(file = stderr(), paste("elpiModulesTable:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/elpiModulesTable.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/elpiModulesTable.RData")
  
  gene_selDF <- as.data.frame(expr_sel$gene_sel)
  rownames(gene_selDF) = gene_selDF[,1]
  gene_selDF = gene_selDF[,-1]
  return(cbind(modules,gene_selDF[modules$feature,]))
})


# temporaImport ----

temporaImport <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaImport started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaImport")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaImport")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaImport", id = "temporaImport", duration = NULL)
  }
  clicked <- input$updatetTemporaParameters
  
  scEx_log <- isolate(scEx_log())
  projections <- isolate(projections())
  tCluster <- isolate(input$temporaCluster)
  tFactor <- isolate(input$temporaFactor)
  tLevels <- isolate(input$temporaLevels)
  
  if (is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), paste("temporaImport:NULL\n"))
    return(NULL)
  }
  
  
  if (!all(c(tCluster,tFactor) %in% colnames(projections))){
    if (DEBUG) cat(file = stderr(), paste("tCluster,tFactor:NULL\n"))
    return(NULL)
  }
  if (!all(tLevels %in% levels(projections[,tFactor]))){
    if (DEBUG) cat(file = stderr(), paste("tLevels:NULL\n"))
    return(NULL)
  }
  if(length(tLevels)<1){
    if (DEBUG) cat(file = stderr(), paste("tLevels:empty\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaImport.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/temporaImport.RData")
  colData(scEx_log) <- S4Vectors::DataFrame(projections[rownames(colData(scEx_log)),])
  # HACK
  # TODO
  # should be a variable, here we are forcing symbol but could be anything.
  rownames(scEx_log) = rowData(scEx_log)[rownames(scEx_log), "symbol"]
  # # tempora hack 
  # suppressMessages(
  #   setMethod("getMD","SingleCellExperiment",
  #             function(x) data.frame(SingleCellExperiment::colData(x)))
  # )
  #in case there are levels where there are no data we have to update this
  colData(scEx_log)[,tFactor] = factor(colData(scEx_log)[,tFactor])
  colData(scEx_log)[,tCluster] = factor(colData(scEx_log)[,tCluster])
  temporaObj <- ImportSeuratObject(seuratobj = scEx_log, 
                                   assayType = "logcounts",
                                   clusters = tCluster,
                                   timepoints = tFactor,
                                   timepoint_order = tLevels,
                                   cluster_labels = levels(projections[,tCluster])
  )
  
  return(temporaObj)
})

# temporaPWProfiles ----
temporaPWProfiles <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaPWProfiles started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaPWProfiles")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaPWProfiles")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaPWProfiles", id = "temporaPWProfiles", duration = NULL)
  }
  BPPARAM=bpparam()
  temporaObj <- temporaImport()
  
  if(is.null(temporaObj)) {
    if (DEBUG) cat(file = stderr(), paste("temporaPWProfiles:NULL\n"))
    return(NULL)
  }
  
  gmt_path = isolate(input$temporaGMTFile)
  min.sz = isolate(input$temporaMinSz)
  max.sz = isolate(input$temporaMaxSz)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaPWProfiles.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/temporaPWProfiles.RData")
  
  if (!file.exists(gmt_path$datapath)) {
    if (DEBUG) cat(file = stderr(), paste("gmt_path:NULL\n"))
    return(NULL)
  }
  
  good = tryCatch ({GSEABase::getGmt(gmt_path$datapath)},
                   error = function(e) {
                     cat(file = stderr(), "GSEABase::getGmt")
                     if (!is.null(getDefaultReactiveDomain())) {
                       showNotification("GSEABase::getGmt not a valid file", 
                                        id = "temporaPWProfiles", type = "error", duration = NULL)
                     }
                     return(NULL)
                   })
  if(is.null(good)) {
    if (DEBUG) cat(file = stderr(), paste("GSEABase::getGmt:NULL\n"))
    return(NULL)
  }
  
  temporaObj <- CalculatePWProfiles(temporaObj, 
                                    gmt_path = gmt_path$datapath,
                                    method="gsva", 
                                    min.sz = min.sz, 
                                    max.sz = max.sz, 
                                    
                                    parallel.sz = BiocParallel::bpnworkers(BPPARAM))
  
  return(temporaObj)
})

# temporaTrajectory ----
temporaTrajectory <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaTrajectory started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaTrajectory")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaTrajectory")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaTrajectory", id = "temporaTrajectory", duration = NULL)
  }
  
  temporaObj <- temporaPWProfiles()
  
  if(is.null(temporaObj)) {
    if (DEBUG) cat(file = stderr(), paste("temporaTrajectory:NULL\n"))
    return(NULL)
  }
  
  n_pcs = isolate(input$temporaNPCs)
  difference_threshold = isolate(input$temporaDiff_thresh)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaTrajectory.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/temporaTrajectory.RData")
  temporaObj = tryCatch ({
    BuildTrajectory(temporaObj, 
                    n_pcs = n_pcs, 
                    difference_threshold = difference_threshold)
  },
  error = function(e) {
    cat(file = stderr(), "tempora BuildTrajectory")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("tempora BuildTrajectory", 
                       id = "temporaBuildTrajectory", type = "error", duration = NULL)
    }
    return(NULL)
  })
  
  
  return(temporaObj)
})


# IdentifyVaryingPWsParallel ----
IdentifyVaryingPWsParallel <- function(object, pval_threshold=0.05){
  require(BiocParallel)
  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }
  if (is.null(object@n.pcs)){
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }
  if (is.null(object@cluster.pathways)){
    stop("CalculatePWProfiles has not been run. See ?Tempora::CalculatePWProfiles for details")
  }
  gsva_bycluster <- object@cluster.pathways
  
  significant_pathways <- c()
  for (i in 1:object@n.pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways <- c(names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5)), significant_pathways)
  }
  
  pca_pathways <- sub("%.*", "", significant_pathways)
  pca_pathways <- gsub("\\s*\\([^\\)]+\\)","",pca_pathways)
  pca_pathways_cleaned <- gsub("[[:punct:]]", "", pca_pathways)
  themes <- pca_pathways_cleaned
  
  if (DEBUG) cat("Fitting GAM models...")
  
  # system.time({
  p_vals <- gams <- list()
  #   for (i in 1:length(themes)){
  #     print(i)
  #     if(length(grep(themes[i], rownames(gsva_bycluster))) == 0) {
  #       p_vals[[i]] <- 1
  #       gams[[i]] <- NA
  #       next
  #     }
  #     if (length(grep(themes[i], rownames(gsva_bycluster))) > 1){
  #       plot_df <- data.frame(cluster=colnames(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=colMeans(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ], na.rm=T))
  #     } else if (length(grep(themes[i], rownames(gsva_bycluster))) == 1){
  #       plot_df <- data.frame(cluster=names(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), value=gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]) }
  #     plot_df$time <- object@cluster.metadata$Cluster_time_score
  #     gams[[i]] <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
  #     temp_anova <- mgcv::anova.gam(gams[[i]])
  #     p_vals[[i]] <- temp_anova$s.pv
  #   }
  # })
  # p_valsOrg = p_vals
  
  func = function(idx, object, themes, gsva_bycluster) {
    require(Tempora)
    require(Matrix)
    require(BiocGenerics)
    # cat(file = stderr(), paste(idx, "\n"))
    grp = BiocGenerics::grep(themes[idx], rownames(gsva_bycluster), ignore.case = T)
    crit = length(grp)
    if(crit == 0) {
      return(list(1, NA))
    }
    if (crit > 1){
      plot_df <- data.frame(
        cluster=colnames(gsva_bycluster[grp, ]), 
        value=Matrix::colMeans(gsva_bycluster[grp, ], 
                       na.rm=T))
    } else if (crit == 1){
      plot_df <- data.frame(cluster=names(gsva_bycluster[grp, ]), 
                            value=gsva_bycluster[grp, ]) 
    } else {
      #should not happen
      return(list(1, NA))
    }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    gams <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
    temp_anova <- mgcv::anova.gam(gams)
    p_vals = temp_anova$s.pv
    list(p_vals, gams)
  }
  
  p_valsNew <- bplapply(1:length(themes), function(x) { func(x, object, themes, gsva_bycluster) })
  for (idx in 1:length(p_valsNew)){
    p_vals[[idx]] = p_valsNew[[idx]][[1]]
    gams[[idx]] = p_valsNew[[idx]][[2]]
  }
  names(p_vals) <- names(gams) <- themes
  
  pval_threshold = pval_threshold
  p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) > 0)]), method = "BH")
  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(names(varying_pathways))]
  
  if (length(varying_pathways)==0){
    cat("No temporally varying pathways detected. Please try running IdentifyVaryingPWs with a more relaxed p-value cutoff.")
    #eventhough the function was not successful return the object because in the vignette
    # this function call sets the original object to what is returned and if it is null
    # you loose all the processing you have done until now.
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("No temporally varying pathways detected.", id = "temporaIdentifyVaryingPWs2", duration = 20)
    }
    return(object)
  } else {
    object@varying.pws <- varying_pathways
    object@gams <- gams
    return(object)
  }
}


# temporaIdentifyVaryingPWs ----
temporaIdentifyVaryingPWs <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaIdentifyVaryingPWs started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaIdentifyVaryingPWs")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaIdentifyVaryingPWs")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaIdentifyVaryingPWs", id = "temporaIdentifyVaryingPWs", duration = NULL)
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetTemporaParameters == 0) {
    return(NULL)
  }
  
  temporaObj <- temporaTrajectory()
  temporaPval_thresh <- input$temporaPval_thresh
  if (is.null(temporaObj) ) {
    if (DEBUG) cat(file = stderr(), paste("temporaIdentifyVaryingPWs:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaIdentifyVaryingPWs.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/temporaIdentifyVaryingPWs.RData")
  
  #Fit GAMs on pathway enrichment profile
  temporaObj <- IdentifyVaryingPWsParallel(object = temporaObj, pval_threshold = temporaPval_thresh)
  
  
  
  return(temporaObj)
})

# temporaPvalModulesTable ----

temporaPvalModulesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "temporaPvalModulesTable started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "temporaPvalModulesTable")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "temporaPvalModulesTable")
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("temporaPvalModulesTable", id = "temporaPvalModulesTable", duration = NULL)
  }
  
  
  # create dependancy on button and return if not pressed once
  if (input$updatetTemporaParameters == 0) {
    return(NULL)
  }
  
  temporaObj <- temporaIdentifyVaryingPWs()
  
  if (is.null(temporaObj) ) {
    if (DEBUG) cat(file = stderr(), paste("temporaPvalModulesTable:NULL\n"))
    return(NULL)
  }
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/temporaPvalModulesTable.RData", list = c(ls()))
  }
  # load(file="~/SCHNAPPsDebug/temporaPvalModulesTable.RData")
  
  outTable = data.frame(goTerm = names(temporaObj@varying.pws))
  outTable$pValues = temporaObj@varying.pws
  
  return(outTable)
  
})


# tempora2DPlotFunc ----
tempora2DPlotFunc <- function(temporaObj, projections, dimX, dimY, dimCol) {
  require(ggnetwork)
  require(ggplot2)
  require(ggrepel)
  require(network)
  space <- projections[, c(dimX, dimY)]
  require(SCORPIUS)
  
  
  # temporaObj@cluster.metadata
  
  data = projections[,c(dimX, dimY, dimCol)]
  if(!all(unlist(lapply(1:2,  FUN = function(x) is.numeric(data[,x]))))) return (NULL)
  
  mean.points <- aggregate(data[, 1:2], list(data[,3]), mean)
  
  # construct a network
  nCl = length(levels(projections[,dimCol]))
  net = matrix(0, nCl, nCl) 
  rownames(net) = levels(projections[,dimCol])
  colnames(net) = levels(projections[,dimCol])
  temporaObj@trajectory$from = as.character( temporaObj@trajectory$from)
  temporaObj@trajectory$to = as.character( temporaObj@trajectory$to)
  apply(temporaObj@trajectory, 1, FUN = function(x){
    net[x[[1]],x[[2]]] <<- 1
    if (! x[[5]] == "unidirectional") {
      net[x[[2]],x[[1]]] <<- 1
    }
  } )
  # traj <- SCORPIUS::infer_trajectory(space)
  n = as.network(net)
  # 
  # gnn = ggnetwork(n, layout = as.matrix(mean.points[,2:3]))
  # # rownames(temporaObj@cluster.metadata) = temporaObj@cluster.metadata$Id
  # 
  # # gnn$vertex.names = temporaObj@cluster.metadata[gnn$vertex.names,"label"]
  # 
  # ggplot(gnn) + geom_edges(aes(x,y,xend = xend,yend=yend),color = "black") + 
  #   geom_nodetext(aes(x,y,label = vertex.names ),
  #                                     fontface = "bold")
  
  dat = ggnetwork(n, layout = as.matrix(mean.points[,2:3]))
  colnames(dat)
  colnames(data) <- c("x", "y", "vertex.names")
  data$xend = NA
  data$yend = NA
  scaleBJ = function(x,y) scale(x, center = min(y), scale = diff(range(y)))
  data$x = scaleBJ(data$x,mean.points[,2])
  data$y = scaleBJ(data$y,mean.points[,3])
  dat = rbind(dat,data)
  p2 = ggplot(dat) +
    geom_point(aes(x, y, col=vertex.names )) +
    geom_edges(data = subset(dat, !is.na(xend)), 
               aes(x, y, xend = xend, yend = yend), 
               color="black",
               alpha = 1, 
               arrow = arrow(length = unit(16, "pt"), type = "closed")) +
    geom_nodes(data = subset(dat, !is.na(xend)), 
               aes(x, y, col=vertex.names), 
               size = 5, color = "white") +
    geom_nodetext(data = subset(dat, !is.na(xend)), 
                  aes(x, y, label = vertex.names), 
                  fontface = "bold") +
    theme_blank() 
  
  p2
}
