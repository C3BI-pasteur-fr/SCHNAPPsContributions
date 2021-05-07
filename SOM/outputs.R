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

# ## Show 'codebook'
# plot(sommap, type="codes", main = "Codes")
# ## Show 'component planes'
# plot(sommap, type = "property", property = sommap$codes[[1]][,1],
#      main = colnames(sommap$codes)[1])
# ## Show 'U-Matrix'
# plot(sommap, type="dist.neighbours")
