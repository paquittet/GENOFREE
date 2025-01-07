


library(ggplot2)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux


source("Fonctions/my_plotabif.R")



ui <- fluidPage(
  fluidRow(
    column(width = 8, class = "well",
           h4("Left plot controls right plot"),
           fluidRow(
             column(width = 6,
                    plotOutput("plot2", height = 300,
                               brush = brushOpts(
                                 id = "plot2_brush",
                                 resetOnNew = TRUE
                               )
                    )
             ),
             column(width = 6,
                    plotOutput("plot3", click = "click_info", height = 300)
             )
           )
    ),
    verbatimTextOutput("info"),
    
    
  )
)

server <- function(input, output) {
  
  # -------------------------------------------------------------------
  # Linked plots (middle and right)
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  
  output$plot2 <- renderPlot({
    ggplot(mtcars, aes(wt, mpg)) +
      geom_point()
  })
  
  output$plot3 <- renderPlot({
    ggplot(mtcars, aes(wt, mpg)) +
      geom_point() +
      coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observe({
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  output$info <- renderText({
    paste0("x=", input$click_info$x, "\ny=", input$click_info$y)
  })
  
  
  
}

shinyApp(ui, server)

