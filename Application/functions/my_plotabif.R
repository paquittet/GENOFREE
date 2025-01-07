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
