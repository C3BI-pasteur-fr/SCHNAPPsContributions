
require(ggplot2)
require(ggalluvial)
# Alluvial plot of two factors
output$alluvial_plot <- renderPlot({
  if (DEBUG) {
    cat(file = stderr(), "alluvial_plot started.\n")
  }
  start.time <- base::Sys.time()
  on.exit({
    printTimeEnd(start.time, "alluvial_plot")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "alluvial_plot")
    }
  })
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("alluvial_plot", id = "alluvial_plot", duration = NULL)
  }
  
  # load reactive data
  projections <- projections()
  alluiv1 <- input$alluiv1
  alluiv2 <- input$alluiv2
  
  # return if nothing to be computed
  if (is.null(projections) ) {
    return(NULL)
  }
  if (alluiv1 == alluiv2 ) {
    return(NULL)
  }
  
  # some debugging messages
  if (DEBUG) cat(file = stderr(), paste("alluvial_plot:\n"))
  # for development and debugging purposes
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/alluvial_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/SCHNAPPsDebug/alluvial_plot.RData")
  
  dat = projections[,c(alluiv1, alluiv2)]
  # dat$cells = rownames(projections)
  gg = ggplot(as.data.frame(dat),
              aes_string( axis1 = alluiv1, axis2 = alluiv2)) +
    geom_alluvium( width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = TRUE) +
    scale_x_discrete(limits = c(alluiv1, alluiv2), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(paste("Alluvial plot of ", alluiv1, "and", alluiv2))
  
  # create and return the plot
  # ggalluvial::(dummyNRow)
  return(gg)
})

# rename projections
observe({
  projections <- projections()
  
  # save(file = "~/SCHNAPPsDebug/alluvial_plot2.RData", list = c(ls(), ls(envir = globalenv())))
  # load(file="~/SCHNAPPsDebug/alluvial_plot2.RData")
  facs = which(lapply(projections, class) == "factor")
  idx=1
  for (ff in facs){
    if(length(levels(projections[,ff]))>100){
      facs = facs[-idx]
    }
    idx = idx + 1
  }
  
  updateSelectInput(session, "alluiv1",
                    choices = c(colnames(projections)[facs])
  )
  updateSelectInput(session, "alluiv2",
                    choices = c(colnames(projections)[facs])
  )
})

