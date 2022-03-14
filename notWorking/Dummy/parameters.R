
#####
# normalization parameters
#####

# myNormalizationChoices ----
# choice for the radio buttion
# the name of the variable will be displayed, the value is the name of the reactive
myNormalizationChoices <- list(
  Dummy_Normalization = "Dummy_Normalization"
)

# myNormalizationParameters ----
# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  Dummy_Normalization = h4("no Parameters implemented")
)

# Dummy_Normalization ----
#' Dummy_Normalization
#' here, we calculate the normalization
#' the result has to be a singleCellExperiment object with the slot logcounts
#' In the example below we devide 
Dummy_Normalization <- reactive({
  # track how much time is spent here
  start.time <- base::Sys.time()
  # remove any notification on exit that we don't want
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "Dummy_Normalization")
  )
  # remove any permanant notification if we rerun reactive
  if (!is.null(getDefaultReactiveDomain()))
    removeNotification(id = "Dummy_NormalizationPerm")
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Dummy_Normalization", id = "Dummy_Normalization", duration = NULL)
  }
  
  
  # just devide by number of cells and scale
  if (DEBUG) cat(file = stderr(), "Dummy_Normalization\n")

  # obviously you only want to use scEx and not scEx_log...
  scEx <- scEx()
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "Dummy_Normalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/Dummy_Normalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/Dummy_Normalization.RData")
  
  
# normalization process
  retVal <- Dummy_NormalizationFunc(scEx = scEx)
# end of normalization. 
  
  addClass("updateNormalization", "green")
  printTimeEnd(start.time, "Dummy_Normalization")
  exportTestValues(Dummy_Normalization = {assays(retVal)[["logcounts"]]})  
  return(retVal)
})


Dummy_NormalizationFunc <- function(scEx, scalingFactor = 10000){
  names(assays(scEx)) <- "logcounts"
  
  return(scEx)
}
