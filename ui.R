#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
# library(plotly)

# Взаимодействие с пользователем
page_fluid(
  title = "Макро-планировщик", 
  theme = bs_theme(version = 5, primary = "#ffa500"), 
  tags$head(tags$link(rel = "stylesheet", href = "styles.css")), 
  
  # Заголовок
  titlePanel(
    uiOutput("title")
  ),
  
  layout_sidebar(
    sidebar = sidebar(
      title = div(
        class = "d-flex gap-2", 
        actionButton(inputId = "btnReset", label = "Сбросить", 
                     class = "btn-outline-secondary btn-sm w-100"), 
        actionButton(inputId = "btnShowAll", label = "Все переменные", 
                     class = "btn-outline-secondary btn-sm w-100")
      ), 
      uiOutput("dials"), 
      width = "30%",
      position = "right",
      open = "desktop"
    ),
    
    navset_card_tab(
      nav_panel(title = "Графики", 
                plotOutput("plots.main", height = "auto"),
                plotOutput("plots.dual", height = "auto")),
      nav_panel(title = "Таблицы", 
                tableOutput("tables.main"), 
                tableOutput("tables.dual")),
      full_screen = T
    ), 
    
    border = F,
    border_radius = F, 
    class = "p-0", 
    position = "right"
  )
)