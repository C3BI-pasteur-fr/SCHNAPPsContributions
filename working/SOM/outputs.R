suppressMessages(require(Rsomoclu))

coE_SOM_dataInput <- callModule(
  cellSelectionModule,
  "coE_SOM_dataInput"
)


# SOM heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSOM",
  coE_heatmapSOMReactive
)

# observer of button Color SOM ----
observe(label = "ob_somParameter", 
        {
          if (DEBUG) cat(file = stderr(), "ob_somParameter\n")
          # browser()
          input$updateSOMParameters
          setRedGreenButtonCurrent(
            vars = list(
              c("coE_geneSOM", input$coE_geneSOM),
              c("coE_dimSOM", input$coE_dimSOM),
              c("coE_SOM_dataInput-Mod_PPGrp", input$'coE_SOM_dataInput-Mod_PPGrp'),
              c("coE_SOM_dataInput-Mod_clusterPP", input$'coE_SOM_dataInput-Mod_clusterPP')
            )
          )
          updateButtonColor(buttonName = "updateSOMParameters", parameters = c(
            "coE_geneSOM", "coE_dimSOM",
            "coE_SOM_dataInput-Mod_PPGrp", "coE_SOM_dataInput-Mod_clusterPP"
          ))
          
        })

observe(label = "somxy",{
  .schnappsEnv$defaultValues[["coE_dimSOMX"]] = input$coE_dimSOMX
  .schnappsEnv$defaultValues[["coE_dimSOMY"]] = input$coE_dimSOMY
})

output$coE_SOMcodebook <- renderPlot({
  sommap = coE_somMapReact()
  if (is.null(sommap)) return(NULL)
  plot(sommap, type="codes", main = "Codes")
})


output$coE_SOMcomponents <- renderPlot({
  sommap = coE_somMapReact()
  if (is.null(sommap)) return(NULL)
  
  plot(sommap, type = "property", property = sommap$codes[[1]][,1],
       main = colnames(sommap$codes)[1])
})
output$coE_SOMuMat <- renderPlot({
  sommap = coE_somMapReact()
  if (is.null(sommap)) return(NULL)
  plot(sommap, type="dist.neighbours")
})

output$coE_somInfo <- renderText({
  
  idxx = input$coE_dimSOMX - 1
  idxy = input$coE_dimSOMY - 1
  res2 = coE_somTrainReact()
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somInfo.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_somInfo.RData")
  
  if (is.null(res2)) {
    return("Genes in neuron\n")
  }
  
  paste("Genes in neuron (", length(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), ")\n",
        paste(names(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), collapse = ", ", sep = ",")
  )
  
})


output$coE_somInfoSymbol <- renderText({
  
  idxx = input$coE_dimSOMX - 1
  idxy = input$coE_dimSOMY - 1
  res2 = coE_somTrainReact()
  scEx_log = scEx_log()
  
  featureData <- rowData(scEx_log)
  geneName = geneName2Index(genesin, featureData)
  
  if (.schnappsEnv$DEBUGSAVE) {
    save(file = "~/SCHNAPPsDebug/coE_somInfo.RData", list = c(ls()))
  }
  # cp = load(file = "~/SCHNAPPsDebug/coE_somInfo.RData")
  
  if (is.null(res2)) {
    return("Genes in neuron\n")
  }
  
  paste("Genes in neuron (", length(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), ")\n",
        paste(featureData[names(which(res2$globalBmus[,1]==idxx & res2$globalBmus[,2]==idxy)), "symbol"], collapse = ", ", sep = ",")
  )
  
})


# ## Show 'codebook'
# plot(sommap, type="codes", main = "Codes")
# ## Show 'component planes'
# plot(sommap, type = "property", property = sommap$codes[[1]][,1],
#      main = colnames(sommap$codes)[1])
# ## Show 'U-Matrix'
# plot(sommap, type="dist.neighbours")
