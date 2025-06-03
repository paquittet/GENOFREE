library(Fragman)
library(ggplot2)
library(gridExtra)

#  Define the folder containing .fsa files 
fsa_folder <- "/Users/medardabiona/Desktop/Genofree/nouvelles_données/89694-M1-090720"

source("/Users/medardabiona/Desktop/Genofree/GENOFREE_fonctions.R")


#  Read all .fsa files in the folder 
my_samples <- storing.inds(folder = fsa_folder)

# #  Define the standard ladder sizes (GS500 LIZ)
# LIZ <- c(75, 100, 139, 150, 160, 200, 300, 350, 400, 450)

############################################################################# #                      Calibration automatique                              # #############################################################################

#  my ladder.info.attach
to.correct <- my_ladder.info.attach(
  stored = my_samples,
  ladder = ladder_used,  # defini dans la foncton de lecture de données standard des fluo + mix à partir d'un fichier
  method = "iter2",
  ladd.init.thresh = NULL,
  channel.ladder = 5,
  draw = TRUE
)

list_data_first_calib <- list.data.covarrubias


#############################################################################
#            Deuxieme calibration des fichiers non calibrés                 #
#############################################################################


# Seuils à tester
thresh_seq <- seq(300, 2000, by = 100)

# Fichiers échoués à la première calibration (corr < 0.999 ou corr ≥ 0.9999)
bad_files <- names(Filter(function(x) {
  is.null(x$corr) || x$corr < 0.999 || x$corr >= 0.9999
}, list_data_first_calib))

# Calibration par fichier
for (file in bad_files) {
  success <- FALSE
  i <- 1
  
  while (!success && i <= length(thresh_seq)) {
    ladd.init.thresh <- thresh_seq[i]
    
    # Calibration du fichier avec un seuil donné
    result <- my_ladder.info.attach(
      stored = my_samples[file],
      ladder = ladder_used,
      method = "iter2",
      ladd.init.thresh = ladd.init.thresh,
      channel.ladder = 3,
      draw = TRUE
    )
    
    # Vérification de la corrélation dans l'environnement global
    if (file %in% names(list.data.covarrubias)) {
      corr_value <- list.data.covarrubias[[file]]$corr
      if (!is.null(corr_value) && corr_value >= 0.999 && corr_value < 0.9999) {
        # Mise à jour uniquement de l'entrée recalibrée
        list_data_first_calib[[file]] <- list.data.covarrubias[[file]]
        success <- TRUE
        
        # Message de confirmation
        message(sprintf("Fichier recalibré : %s avec seuil %d (corr = %.4f)", file, ladd.init.thresh, corr_value))
      } else {
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  if (!success) {
    message(sprintf("Calibration échouée pour %s après tous les seuils testés", file))
  }
}








