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


