library(shiny)
library(seqinr)
library(ggplot2)



# FONCTIONS --------------------------------------------------------
my_plotabif <- function(abifdata = data,
                        n_ticks,  # add control for tick number
                        chanel = 1,
                        tmin = 1/tscale,
                        tmax = abifdata$Data[["SCAN.1"]]/tscale,
                        tscale = 1000,
                        yscale = 1000,
                        xlab = paste("Time", tscale, sep = "/"),
                        ylab = paste("RFU", yscale, sep = "/"),
                        irange = (tmin * tscale):(tmax * tscale),
                        x = irange/tscale,
                        xlim = NA,
                        chanel.names = c(1:4, 105),
                        DATA = paste("DATA", chanel.names[chanel], sep = "."),
                        y = abifdata$Data[[DATA]][irange]/yscale,
                        ylim = c(0, max(y)),  # change min(y) to 0
                        main = paste0("Fluochrome : ", abifdata$Data[[paste("DyeN", chanel, sep = ".")]]),
                        calibr = NULL,
                        ...) {
  
  # Change color according to fluochrome
  if (chanel == 1) {
    color = "#1C86EE"
  } else if (chanel == 2) {
    color = "#00CD00"
  } else if (chanel == 3) {
    color = "#EEC900"
  } else if (chanel == 4) {
    color = "#EE2C2C"
  } else {
    color = "black"
  }
  
  # Apply calibration if provided (ensuring calibration only affects x)
  x <- calibr(irange)  # Ensure calibr is applied consistently
  
  # Set xlim if NA
  if (any(is.na(xlim))) {
    xlim <- c(0, max(x))  # Adjust xlim based on calibrated x values
  }
  
  # Create data frame for ggplot
  df_gg <- data.frame(x = x, y = y)
  
  # PLOT
  plot_gg <- ggplot(df_gg, mapping = aes(x = x, y = y)) +
    geom_line(color = color) +
    
    # Optionally set tick marks if required (uncomment if needed)
    # scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = as.numeric(n_ticks))) +
    
    # Use a default ggplot theme
    theme_minimal() +
    
    # Label xlab, ylab and title
    labs(x = xlab,
         y = ylab,
         title = main) +
    
    # Set y and x limits (ensure scaling works with calibration)
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    
    theme(
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      plot.title = element_text(size = 18, margin = margin(b = 6)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.y = element_text(size = 18, margin = margin(0, 9, 0, 0), face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.text.x = element_text(margin = margin(t = 6), angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(margin = margin(r = 6)),
      text = element_text(size = 18)
    )
  
  return(plot_gg)  # Return the plot directly
}

# USER INTERFACE ----------------------------------------------------------
# Define the UI of the app
ui <- fluidPage(
  titlePanel(h1("VISUALISATION DES CHROMATOGRAMMES")),
  
  sidebarLayout(
    sidebarPanel(
      # IMPORT BUTTON
      fileInput(inputId = "fsa_file",  
                label = "Fichier .fsa",
                multiple = FALSE,
                accept = ".fsa",
                placeholder = "Sélectionnez le fichier .fsa",
                buttonLabel = "Parcourir"),
      
      # RECALIBRAGE BUTTON (tmin)
      numericInput(
        inputId = "tmin",
        label = h3("Recalibrage"),
        value = 1.5,
        min = 1,
        max = 4,
        step = 0.05,
      ),
      
      # FLUOCHROME CHOICE
      radioButtons(inputId = "fluo_choice",
                   label = h3("Fluochrome"),
                   choices = list("FAM" = 1, "VIC" = 2, "NED" = 3, "PET" = 4, "CONTROL" = 5),
                   selected = 1),
      
      # XLIM OF CHROMMATOGRAMME
      sliderInput(
        inputId = "xlim",
        label = h3("Intervalle de l'axe des x"),
        value = c(0, 300),
        min = 0,
        max = 300,
        step = 10,
      ),
      
      # YLIM OF CHROMMATOGRAMME
      sliderInput(
        inputId = "ylim",
        label = h3("Intervalle de l'axe des y"),
        value = c(0, 30),
        min = 0,
        max = 50,
        step = 5,
      ),
      
      # NUMBER OF TICKS N_TICKS
      selectInput(inputId = "nticks",
                  multiple = F,
                  label = h3("Précision de l'axe des x"), 
                  choices = list("0.5" = 0.5, "1" = 1, "5" = 5, "10" = 10, "50" = 50), 
                  selected = 10)
      
    ),
    
    mainPanel(
      fluidRow(
        column(6, plotOutput("plot_fsa", click = "plot_click", brush = brushOpts(id = "plot_fsa_brush", resetOnNew = T))),
        column(6, plotOutput("plot_fsa_zoom", click = "plot_zoom_click"))
      ),
      verbatimTextOutput("info"),
      verbatimTextOutput("info2"),
      
      plotOutput("plot_peak")
    )
  )
)









# SERVER -----------------------------------------------------------------
server <- function(input, output, session) {
  
  ###########################################################################
  # REACTIVE VALUES  ------------------------------------------------------------------
  ###########################################################################
  
  # VARIABLES OUTSIDE FUNCTION (pour éviter de les préciser à chaque fois)
  # Use a reactiveVal to store the .fsa data
  fsa_data_reactive <- reactiveVal(NULL)
  
  # When a file is uploaded, read it and store the data in fsa_data
  observeEvent(input$fsa_file, {
    req(input$fsa_file)
    fsa_data_reactive(read.abif(input$fsa_file$datapath))
  })
  
  
  
  
  ###########################################################################
  # OUTPUT ------------------------------------------------------------------
  ###########################################################################
  
  # 1. PLOT PEAKS
  output$plot_peak <- renderPlot({
    req(fsa_data_reactive())  # Ensure the file is loaded
    data <- fsa_data_reactive()  # Get the data from the reactive value
    
    # Peak location calculation (channel 5 is the standard)
    maxis <- peakabif(
      data,
      chanel = 5,
      npeak = 14,  # Standard with 14 peaks
      tmin = input$tmin
    )
    
    # Calibration data
    data(gs500liz)
    lizok <- gs500liz$liz
    lizok[!gs500liz$mask1 | !gs500liz$mask2] <- NA
    lizok <- lizok[-c(1,2)]
    y <- lizok[!is.na(lizok)]
    x <- maxis$maxis[!is.na(lizok)]
    
    # Calibration function
    calibrage_function <- splinefun(x, y)  # local evaluated the expression here in shiny
  })
  
  
  
  
  
  # 2. PLOT CHROMATOGRAMME
  output$plot_fsa <-   renderPlot({
    # PLOT CHROMATOGRAMME
    req(fsa_data_reactive())  # Ensure the file is loaded
    data <- fsa_data_reactive()  # Get the data from the reactive value
    
    # CALIBRATION FUNCTION depending on Tmin value
    maxis <- peakabif(
      data,
      chanel = 5,
      npeak = 14,  # Standard with 14 peaks
      tmin = input$tmin,
      fig = F
    )
    
    # Calibration data
    data(gs500liz)
    lizok <- gs500liz$liz
    lizok[!gs500liz$mask1 | !gs500liz$mask2] <- NA
    lizok <- lizok[-c(1,2)]
    y <- lizok[!is.na(lizok)]
    x <- maxis$maxis[!is.na(lizok)]
    
    # Calibration function
    calibrage_function <- splinefun(x, y)
    
    
    # PLOT NO ZOOM
    my_plotabif(
      abifdata = data,
      chanel = as.numeric(input$fluo_choice),
      calibr = calibrage_function,
      xlim = c(as.numeric(input$xlim)),
      ylim = c(as.numeric(input$ylim)),
      n_ticks = as.numeric(as.character(input$nticks))
    )  # local evaluated the expression here in shiny
  })
  
  
  
  
  
  
  
  
  

  
  # 3. PLOT ZOOM
  # REACTIVE ZOOM RANGE
  ranges <- reactiveValues(x = NULL, y = NULL)  # for zooming
  
  output$plot_fsa_zoom <- renderPlot({
    # PRELIMINARY CHECKS
    req(fsa_data_reactive())  # Ensure the file is loaded
    data <- fsa_data_reactive()  # Get the data from the reactive value
    req(ranges$x, ranges$y)  # Explicitly react to ranges, fait disparaitre le plot zoomé si aucun brush
    
    
    # CALIBRATION FUNCTION depending on Tmin value
    maxis <- peakabif(
      data,
      chanel = 5,
      npeak = 14,  # Standard with 14 peaks
      tmin = input$tmin,
      fig = F
    )
    
    # Calibration data
    data(gs500liz)
    lizok <- gs500liz$liz
    lizok[!gs500liz$mask1 | !gs500liz$mask2] <- NA
    lizok <- lizok[-c(1,2)]
    y <- lizok[!is.na(lizok)]
    x <- maxis$maxis[!is.na(lizok)]
    
    # Calibration function
    calibrage_function <- splinefun(x, y)
    
    # Plotting the chromatogram with zoom applied
    plot_gg <- 
      my_plotabif(
        abifdata = data,
        chanel = as.numeric(input$fluo_choice),
        calibr = calibrage_function,
        xlim = c(as.numeric(input$xlim)),
        ylim = c(as.numeric(input$ylim))
      )
    
    suppressWarnings(print(
      plot_gg +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = T)
    ))
  })
  
  # ENABLE ZOOM
  observe({
    brush <- input$plot_fsa_brush  
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  
  
  # 4. COORDINATES OF THE MOUSE
  output$info <- renderText({
    paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
  })
  
  output$info2 <- renderText({
    paste0("x=", input$plot_zoom_click$x, "\ny=", input$plot_zoom_click$y)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
