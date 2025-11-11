library(CelliD)
# max expressed genes ----
callModule(
  tableSelectionServer,
  "cellID_cellTable",
  cellID_cellTableReact
)

callModule(
  tableSelectionServer,
  "cellID_geneTable",
  cellID_geneTableReact
)


callModule(
  tableSelectionServer,
  "cellID_stdevTable",
  cellID_stdevTableReact
)

# change color of button
observe(label = "updateCellID", {
  if (DEBUG) cat(file = stderr(), "observe updateCellID\n")
  out <- cellIdReactive()
  if (is.null(out)) {
    .schnappsEnv$calculated_gQC_tsneDim <- "NA"
  }
  input$updateCellID
  
  setRedGreenButtonCurrent(
    vars = list(
      c("cellID_Method", input$cellID_Method),
      c("cellID_nmcs", input$cellID_nmcs)
    )
  )
  updateButtonColor(buttonName = "updateCellID", parameters = c(
    "cellID_Method", "cellID_nmcs"
  ))
})

# change projections to use 

observe(label = "cellID_gBy",{
  projFactors <- projFactors()
  if (is.null(projections)) {
    return(NULL)
  }
  updateSelectInput(
    session,
    "cellID_gBy",
    choices = projFactors
  )
})

callModule(
  tableSelectionServer,
  "cellID_GeneSetTable",
  cellIDGeneSetReact
)


output$cellID_stdev <- renderPlot({
  require(ggplot2)
  if (DEBUG) cat(file = stderr(), "cellID_stdev started.\n")
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "cellID_stdev")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellID_stdev")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("cellID_stdev", id = "cellID_stdev", duration = NULL)
  }
  mca = cellIdReactive()
  if (is.null(mca)) {
    return(NULL)
  }
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/cellID_stdev.RData", list = c(ls()))
  }
  # load(file = "~/SCHNAPPsDebug/cellID_stdev.RData")
  
  # h2("Variances of PCs")
  
  
  
  df <- data.frame(stdev = attr(reducedDim(mca, "MCA"), "stdev"))
  df$component <- as.numeric(rownames(df))
  retVal <- ggplot(data = df,aes(x=component, y=stdev)) + geom_bar(stat = "identity")  
  exportTestValues(cellID_stdev = {
    retVal
  })
  return(retVal)
  # barplot(pca$var_pcs, main = "Variance captured by first PCs")
})


