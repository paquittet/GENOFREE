library(shiny)
library(seqinr)
library(ggplot2)
library(shinythemes)
library(tidyverse)
library(shinyjs)
library(openxlsx)
library(rhandsontable)
library(bslib)
library(Fragman)



# Define UI for app that draws a histogram ----
ui <- page_sidebar(
  # App title ----
  title = "Hello Shiny!",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    # Input: Slider for the number of bins ----
    sliderInput(
      inputId = "bins",
      label = "Number of bins:",
      min = 1,
      max = 50,
      value = 30
    )
  ),
  # Output: Histogram ----
  plotOutput(outputId = "distPlot")
)









server <- function(input, output) {
   }











# Run the application 
shinyApp(ui = ui, server = server)