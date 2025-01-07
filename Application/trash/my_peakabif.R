abifdata = data
tscale = 1000
yscale = 1000
chanel = 5
npeak = 14
thres = 400 / yscale
fig = TRUE
chanel.names = c(1:4, 105)
DATA = paste("DATA", chanel.names[chanel], sep = ".")
tmin = 1
tmax = abifdata$Data[["SCAN.1"]] / tscale
irange = (tmin * tscale):(tmax * tscale)
y = abifdata$Data[[DATA]][irange] / yscale
method = "monoH.FC"
maxrfu = 1000


# function (abifdata,
#           chanel,
#           npeak,
#           thres = 400 / yscale,
#           fig = TRUE,
#           chanel.names = c(1:4, 105),
#           DATA = paste("DATA", chanel.names[chanel], sep = "."),
#           tmin = 1 / tscale,
#           tmax = abifdata$Data[["SCAN.1"]] / tscale,
#           tscale = 1000,
#           yscale = 1000,
#           irange = (tmin * tscale):(tmax *
#                                       tscale),
#           y = abifdata$Data[[DATA]][irange] / yscale,
#           method = "monoH.FC",
#           maxrfu = 1000,
#           ...)
# 
# {

tmin = 1
tmax = abifdata$Data[["SCAN.1"]] / tscale
irange = (tmin * tscale):(tmax * tscale)
y = abifdata$Data[[DATA]][irange] / yscale


# Mettre à zéro tous les points dont l'intensité est inférieure au seuil (thres)
y[y < thres] <- 0

# Initialiser les vecteurs pour stocker les résultats des différents paramètres des pics
heights <- surfaces <- maxis <- starts <- stops <- numeric(npeak)

# Variable pour savoir si nous sommes dans une zone de bruit (initialisé à TRUE)
innoise <- TRUE

# Compteur d'index pour suivre le numéro du pic (initialisé à 1)
pkidx <- 1

# Boucle pour parcourir toutes les valeurs de y et détecter les pics
for (i in 1:length(y)) {
  if (y[i] > 0) {  # Si l'intensité est positive (indiquant un pic)
    if (innoise) {  # Si on était dans une zone de bruit
      starts[pkidx] <- i  # Enregistrer la position de début du pic
      innoise <- FALSE  # Passer à l'état de "pic" (on n'est plus dans le bruit)
    }
  }
  else {  # Si l'intensité est égale à zéro (ou négative, ce qui est peu probable)
    if (!innoise) {  # Si on sort d'une zone de pic
      stops[pkidx] <- i - 1  # Enregistrer la position de fin du pic (exclure le zéro)
      innoise <- TRUE  # Passer à l'état de "bruit"
      pkidx <- pkidx + 1  # Passer au prochain pic
    }
  }
}

# Si 'fig' est TRUE, afficher plusieurs graphiques dans une grille 4x4
if (fig) 
  graphics::par(mfrow = c(4, 4), mar = c(2, 2, 0, 0) + 0.2, oma = c(0, 0, 2, 0))

# Boucle pour analyser chaque pic détecté
for (i in 1:npeak) {
  x <- starts[i]:stops[i]  # Extraire les indices du pic actuel, de start à stop
  if (length(x) <= 2) {  # Si le pic est trop petit (moins de 3 points), le sauter
    maxis[i] <- NA
    warning("Not all requested peaks were assigned")  # Avertir qu'un pic a été ignoré
    next  # Passer à l'itération suivante de la boucle
  }
  
  # Créer une fonction spline pour approximer la courbe du pic
  spfun <- stats::splinefun(x, y[x], method = method)
  
  # Trouver le maximum du pic en optimisant la fonction spline
  maxis[i] <- stats::optimize(spfun, interval = range(x), maximum = TRUE)$maximum
  
  # Calculer la hauteur du pic à la position du maximum
  heights[i] <- spfun(maxis[i])
  
  # Calculer l'aire sous la courbe (surface) du pic en intégrant la fonction spline
  surfaces[i] <- stats::integrate(spfun, starts[i], stops[i])$value
  
  # Si 'fig' est TRUE, afficher un graphique pour ce pic
  if (fig) {
    xx <- (x - 1)/tscale + tmin  # Ajuster les indices pour l'échelle du temps
    graphics::plot(xx, y[x], type = "p", las = 1, ylim = range(y))  # Afficher les points du pic
    graphics::abline(h = thres, col = "red")  # Ajouter une ligne rouge représentant le seuil
    graphics::lines(xx, spfun(x), col = "blue")  # Tracer la courbe spline du pic en bleu
    graphics::abline(v = (maxis[i] - 1)/tscale + tmin, col = "grey")  # Ajouter une ligne verticale pour le maximum du pic
  }
}

# Estimer la ligne de base (bruit) à partir des données ABIF
baseline <- baselineabif(abifdata$Data[[DATA]][irange], maxrfu = maxrfu)

# Ajuster la ligne de base en fonction de l'échelle de l'axe des ordonnées (yscale)
baseline <- baseline/yscale

# Si 'fig' est TRUE, afficher un texte d'information sur les paramètres utilisés
if (fig) 
  graphics::mtext(paste(deparse(substitute(abifdata)), 
                        ",", DATA, ", tmin =", tmin, ", tmax =", tmax, ", thres =", 
                        thres, ", npeak =", npeak, ", yscale = ", yscale), 
                  side = 3, outer = TRUE)

# Retourner les résultats sous forme de liste, avec les pics analysés et la ligne de base ajustée
invisible(list(maxis = (maxis - 1) + tmin * tscale, 
               heights = yscale * heights, 
               surfaces = yscale * surfaces, 
               baseline = baseline))
