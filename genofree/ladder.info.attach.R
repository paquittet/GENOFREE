library(Fragman)
library(ggplot2)
library(gridExtra)

#  Define the folder containing .fsa files 
fsa_folder <- "~/Desktop/Genofree/genofree/fichier .fsa test/Re Fichier fsa test/fsa"

#  Read all .fsa files in the folder 
my.samples <- storing.inds(folder = fsa_folder)

#  Define the standard ladder sizes (GS500 LIZ)
LIZ <- c(75, 100, 139, 150, 160, 200, 300, 350, 400, 450)

calib_env <- new.env()


#  ladder.info.attach
to.correct <- ladder.info.attach(
  stored = my.samples,
  ladder = LIZ,
  method = "iter2",
  ladd.init.thresh = 1500,
  channel.ladder = 5,
  draw = TRUE,
  prog = TRUE,
  env = calib_env
)
 # Extract the calibration results
calib_results <- calib_env$list.data.covarrubias

sample_name <- "A02_1890_016.fsa" # just 1 .fsa file
calib_data <- calib_results[[sample_name]]

# Extract the ladder peak positions, heights, and weights

pos <- calib_data$pos # Time
hei <- calib_data$hei # Peak
wei <- calib_data$wei # ladder size(bp)

# plot 1 (chromatogram with ladder peaks)

chromatogram_df <- data.frame(
  Time = pos / 1000,        # Scale down for display
  RFU = hei / 1000          # Scale down for display
)

plot1 <- ggplot(chromatogram_df, aes(x = Time, y = RFU)) +
  geom_segment(aes(xend = Time, yend = 0), color = "black", linewidth = 1.2) +
  labs(title = "Ladder Peaks Chromatogram",
       x = "Time / 1000", y = "RFU / 1000") +
  theme_minimal()

# Plot 2: ladder size vs migration time
calib_df <- data.frame(
  Size = calib_data$wei,
  Time = calib_data$pos
)

plot2 <- ggplot(calib_df, aes(x = Size, y = Time)) +
  geom_point(size = 3, color = "black", alpha = .8) +
  geom_line(color = "black", alpha = .8) +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed") +
  labs(
    title = "Calibration curve (Expected LIZ sizes vs. Migration time)",
    x = "Standard size (bp)",  # LIZ
    y = "Migration time"       # pos
  ) +
  theme_minimal()

plot2 <- ggplot(calib_df, aes(x = Size, y = Time)) +
  geom_point(size = 3, color = "black") +
  geom_line(color = "black") +
  geom_text(aes(label = Size), vjust = -1, size = 3.5, color = "#7CCD7C", alpha = .8) +  
  geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed", alpha = .8) +
  labs(
    title = "Calibration  withLIZ",
    x = "Standard size (bp)",  # LIZ
    y = "Migration time"
  ) +
  theme_minimal()

# display side by side 

grid.arrange(plot1, plot2, ncol = 2)