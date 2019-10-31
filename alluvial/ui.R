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
  tabItem("alluvialTab",
          box(
            title = "alluvial plot", solidHeader = TRUE, width = 12, status = 'primary', 
            fluidRow(
              column(width = 6, 
                     selectInput("alluiv1", "1st axsis", choices = c("notyet"), selected = "notyet")),
              column(width = 6, 
                     selectInput("alluiv2", "2nd axsis", choices = c("notyet"), selected = "notyet"))
            ),
            fluidRow(
              column(width = 12, 
                     plotOutput("alluvial_plot") # %>% withSpinner()
              )
            )
          )
          
  )
)

