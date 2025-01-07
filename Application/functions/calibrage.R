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
calibrage_function <- splinefun(x, y)
