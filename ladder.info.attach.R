################################################################################
#                                                                              #
#         CALIBRATION AUTOMATIQUE DES FICHIERS .FSA                            #
#                                                                              #
################################################################################
#' *CONTEXTE GÉNÉRAL* :
#' Ce script permet d'effectuer une calibration automatique des fichiers 
#' chromatogrammes (.fsa) issus de génotypages microsatellites. 
#' Il repose sur les fonctions du package **Fragman** ainsi que sur des fonctions 
#' personnalisées définies dans le fichier `GENOFREE_fonctions.R`.
#'
#' *OBJECTIFS* :
#' - Charger et lire les fichiers .fsa à partir d’un dossier donné
#' - Effectuer une première calibration avec un seuil automatique
#' - Identifier les fichiers mal calibrés
#' - Appliquer une recalibration itérative uniquement sur ceux-ci
#'
#' *DÉPENDANCES* :
#' - Fragman : pour la lecture et la calibration des fichiers .fsa
#' - ggplot2, gridExtra : pour les visualisations
#'
#' *ENTRÉES* :
#' - `fsa_folder` : chemin du dossier contenant les fichiers .fsa
#' - `GENOFREE_fonctions.R` : contient les versions modifiées de fonctions Fragman

library(Fragman)
library(ggplot2)
library(gridExtra)
library(commentr)

#  Define the folder containing .fsa files 

# fsa_folder <- "/data"

fsa_folder <- "/Users/medardabiona/Desktop/Genofree/Mix1"

# source("GENOFREE_fonctions.R")

source("/Users/medardabiona/Desktop/Genofree/GENOFREE_fonctions.R")

#  Read all .fsa files in the folder 
my_samples <- storing.inds(folder = fsa_folder)


################################################################################
#                                                                              #
#       PREMIÈRE CALIBRATION DES FICHIERS .FSA AVEC MY LADDER.INFO.ATTACH      #
#                                                                              #
################################################################################
#' *CONTEXTE* :
#' Cette étape applique une première calibration à l’ensemble des fichiers .fsa 
#' en utilisant la fonction `my_ladder.info.attach()` (version modifiée de 
#' `ladder.info.attach()` du package Fragman).
#' ensuite stock le resultat issue de cette calibration dans l'env 
#' list.data.covarrubias
#' 
#' *REMARQUE* :
#' `ladd.init.thresh = NULL` laisse la fonction choisir automatiquement un seuil 
#' initial de détection. Cette calibration peut échouer sur certains fichiers 
#' qui seront ensuite traités via une procédure de recalibration.

to.correct <- my_ladder.info.attach(
  stored = my_samples,
  ladder = ladder_used,  #' defini dans la foncton de lecture de données standard 
  #' des fluo + mix à partir d'un fichier
  method = "iter2",
  ladd.init.thresh = NULL,
  channel.ladder = 5,
  draw = TRUE
)

# Résultat initial de calibration
list_data_first_calib <- list.data.covarrubias





################################################################################
#                                                                              #
#              RECALIBRATION AUTOMATIQUE DES FICHIERS MAL CALIBRÉS             #
#                                                                              #
################################################################################
#' *CONTEXTE* :
#' Cette étape intervient après la première calibration, pour les fichiers 
#' dont la corrélation est jugée insuffisante (< 0.999 ou > 0.9999).
#' Elle repose sur des itérations avec différents seuils `ladd.init.thresh`.
#'
#' *OBJECTIF* :
#' - Identifier les fichiers échoués à la première calibration
#' - Tester automatiquement plusieurs seuils de calibration de 300 à 2000
#' - Valider la calibration si la corrélation obtenue est dans l’intervalle
#'  [0.999 ; 0.9999]
#'
#' *TRAITEMENT* :
#' - La calibration est tentée seuil par seuil jusqu’à réussite
#' - Les calibrations réussies sont mises à jour dans `list.data.covarrubias`
#'
#' *REMARQUE* :
#' Seuls les fichiers mal calibrés sont recalibrés ; ceux avec une corrélation 
#' déjà comprise entre 0.999 et 0.9999 ne sont pas retraités.


# Seuils à tester
# Paramètres de seuils
corr_min <- 0.999
corr_max <- 0.99994
thresh_seq <- seq(300, 2000, by = 100)

# Fichiers échoués à la première calibration (corr < corr_min ou corr > corr_max)
bad_files <- names(Filter(function(x) {
  is.null(x$corr) || x$corr < corr_min || x$corr > corr_max
}, list.data.covarrubias))

################################################################################
#                                                                              #
# La boucle for parcourt tous les fichiers mal calibrés (bad_files) et tente,  #
# pour chacun, de les recalibrer automatiquement en testant une série de       #
# seuils (thresh_seq).                                                         #
#                                                                              #
################################################################################

for (file in bad_files) { 
  success <- FALSE  # Indique si la calibration a réussi pour ce fichier
  i <- 1.     # Index de seuil de calibration à tester dans thresh_seq
  
  # teste plusieurs seuils jusqu'à obtenir une bonne calibration ou épuiser la liste
  while (!success && i <= length(thresh_seq)) {
    ladd.init.thresh <- thresh_seq[i]
    
    # Calibration du fichier avec un seuil donné
    result <- my_ladder.info.attach(
      stored = my_samples[file],
      ladder = ladder_used,
      method = "iter2",
      ladd.init.thresh = ladd.init.thresh,
      channel.ladder = 5,
      draw = TRUE
    )
    
    # Vérification de la corrélation dans l'environnement global
    if (file %in% names(list.data.covarrubias)) {
      corr_value <- list.data.covarrubias[[file]]$corr
      if (!is.null(corr_value) && corr_value >= corr_min && corr_value <= corr_max) {
        # Mise à jour uniquement de l'entrée recalibrée
        list_data_first_calib[[file]] <- list.data.covarrubias[[file]]
        success <- TRUE
        
        # Message de confirmation
        message(sprintf("Fichier recalibré : %s avec un threshold de %d (corr = %.4f)", file, ladd.init.thresh, corr_value))
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




