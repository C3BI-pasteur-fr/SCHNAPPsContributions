# menu list
# defines the main entry
menuList <- list(
  menuItem("Alluvial Plots",
    tabName = "alluvialMenuList", startExpanded = FALSE,
    # id = "dummyID",
    menuSubItem("alluvialTab", tabName = "alluvialTab")
  )
)


# tabs with the actual content
tabList <- list(
  crHeatMapTab = tabItem(
    "alluvialTab",
    fluidRow(
      column(4, offset = 1,
             selectInput("alluiv1", "1st axsis", choices = c("notyet"), selected = "notyet")),
    ),
    fluidRow(
      column(4, offset = 0,
             selectInput("alluiv2", "2nd axsis", choices = c("notyet"), selected = "notyet")),
    ),
    tags$h3("alluvial plot"),
    fluidRow(
      plotOutput("alluvial_plot") # %>% withSpinner()
    )
  )
)

