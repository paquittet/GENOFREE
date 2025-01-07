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
