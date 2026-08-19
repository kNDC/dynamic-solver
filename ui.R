## Пользовательская составляющая
# https://shiny.posit.co/

library(shiny)
library(bslib)
# library(plotly)

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
      nav_panel(title = "Графики", div(
                  plotOutput("plots_main", height = "auto"),
                  plotOutput("plots_dual", height = "auto"))),
      nav_panel(title = "Таблицы", 
                tableOutput("tables_main"), 
                tableOutput("tables_dual")),
      full_screen = TRUE
    ), 
    
    border = FALSE,
    border_radius = FALSE, 
    class = "p-0", 
    position = "right"
  )
)