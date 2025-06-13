################################################################################
#                                                                              #
#                         FONCTIONS FRAGMAN MODIFIEES                          #
#                                                                              #
################################################################################
#' *CONTEXTE* : Fragman ne permet pas d'exploiter directement les sorties de 
#' fonctions, donc on modifie les fonctions du package de façon à contrôler les
#' élements d'output des fonctions ladder.info.attach(), finder.ladder()...










################################################################################
#                                                                              #
#                                PLOT_CALIBR()                                 #
#                                                                              #
################################################################################
#' *DESCRIPTION*
#' Cette fonction plot le chromatogrammes de la channel 5 (le standard) et les 
#' pics détectés automatiquement en rouge. C'est la partie graphique de la fonction
#' find.ladder() modifiée. Elle est donc intégrée dans cette dernière, elle-même
#' intégrée dans la fonction ladder.info.attach().

#' *ARGUMENT* 
#' @param x - les données de fluroescence du chanel précisé dans ladder.info.attach
#' @param roxy3 - contient les coordonnées des pics détectés le R2
#' @param limi - limite de l'axe des ordonnées, calculés dans ladder.info.attach()
# plot_calibr <- function(x = x, roxy3 = roxy3, limi = limi){
#   data_gg <- data.frame(rfu = as.numeric(x), index = 1:length(x))
#   points_df <- data.frame(
#     pos = roxy3$pos,
#     hei = roxy3$hei,
#     type = "Peaks selected"  # Cela servira à nommer la légende
#   )
#   corr_text <- paste("Correlation:", round(roxy3$corr, digits = 4))
#   
#   plot_temp <- ggplot(data = data_gg, aes(y = rfu, x = index)) +
#     geom_line(lwd = 0.9, col = transp("black", 0.8)) +
#     ylim(c(0, (limi[3] + 1000))) +
#     ylab("RFU") +
#     xlab("") +
#     ggtitle(attributes(x)$mycomm) +
#     theme_classic(base_size = 19) +
#     scale_x_continuous(
#       breaks = roxy3$pos,
#       labels = roxy3$wei
#     ) +
#     geom_point(
#       data = points_df,
#       aes(x = pos, y = hei),
#       size = 5,
#       shape = 19
#     ) +
#     geom_point(
#       data = points_df,
#       aes(x = pos, y = hei),
#       size = 3,
#       shape = 19,
#       col = "red"
#     ) +
#     scale_color_manual(
#       values = c("Peaks selected" = "red"),
#       name = NULL
#     ) +
#     # Ajouter une annotation en haut à gauche pour la "légende" de la corrélation
#     annotate("label",
#              x = min(data_gg$index),
#              y = limi[3] + 950,
#              label = corr_text,
#              fontface = "bold",
#              hjust = 0,
#              color = "red",
#              size = 7,
#              label.size = 0.3,      # épaisseur du cadre
#              label.r = unit(0.15, "lines"),  # arrondi des coins
#              fill = "white")  +      # couleur de fond du label
#     theme(
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
#     )
#   
#   return(plot_temp)
# }

plot_calibr <- function(x = x, roxy3 = roxy3, limi = limi) {
  # 1. Données du signal
  data_gg <- data.frame(rfu = as.numeric(x), index = 1:length(x))
  
  # 2. Points détectés (pics sélectionnés)
  points_df <- data.frame(pos = roxy3$pos, hei = roxy3$hei)
  
  # 3. Texte de corrélation
  # corr_text <- paste("Correlation:", round(roxy3$corr, 4))
  
  corr_value <- if (is.null(roxy3$corr) || is.na(roxy3$corr) || roxy3$corr < 0.999 || roxy3$corr > 0.99994) {
    "N/A"
  } else {
    sprintf("%.4f", roxy3$corr)
  }
  corr_text <- paste("Correlation:", corr_value)
  
  # 4. Définir une échelle Y automatique (pour éviter d’écraser les pics)
  ylim_top <- max(c(data_gg$rfu, roxy3$hei), na.rm = TRUE) * 1.15
  
  # 5. Création du graphique
  plot_temp <- ggplot(data = data_gg, aes(x = index, y = rfu)) +
    
    # Courbe du chromatogramme
    geom_line(linewidth = 0.8, color = "black") +
    
    # Points rouges (pics sélectionnés)
    geom_point(data = points_df,
               aes(x = pos, y = hei),
               size = 3.2,
               shape = 19,
               color = "#8B0000") +
    
    # Affichage dynamique de la corrélation
    annotate("label",
             x = min(data_gg$index) + 150,
             y = ylim_top - 0.04 * ylim_top,
             label = corr_text,
             fontface = "bold",
             hjust = 0,
             color = "#8B0000",
             size = 5,
             label.size = 0.25,
             label.r = unit(0.2, "lines"),
             fill = "white") +
    
    # Axes
    ylim(0, ylim_top) +
    xlab("Temps") +
    ylab("RFU") +
    
    # Titre dynamique (nom du fichier)
    ggtitle(paste("Chromatogramme", attributes(x)$mycomm)) +
    
    # Style graphique
    theme_classic(base_size = 18) +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_line(color = "grey95", size = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(color = "black", size = 0.6),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 18, margin = margin(b = 10))
    )
  
  return(plot_temp)
}



################################################################################
#                                                                              #
#                               MY_FIND.LADDER()                               #
#                                                                              #
################################################################################
#' *DESCRIPTION*
#' Fonction modifiée de la fonction Fragman::find.ladder(), qui inclut plot_calibr()
#' à la place de la partie graphique originelle. Tout le reste de la fonction est 
#' identique à l'originale. Les modifications sont entourés de blocs de 
#' commentaires "MODIFIED"

#' *ARGUMENT* 
#' @param : identique à la fonction find.ladder() originale

my_find.ladder <- 
  function (x,
          ladder,
          draw = TRUE,
          dev = 50,
          warn = TRUE,
          init.thresh = NULL,
          sep.index = 8,
          method = NULL,
          reducing = NULL,
          who = "sample",
          attempt = 10,
          cex.title = 0.8)

{
  hohoho <- big.peaks.col(x[1:length(x)], 100)
  if (is.null(method)) {
    tototo <- big.peaks.col(x[1:length(x)], median(hohoho$hei) * 
                              2)
    close <- (length(x) - tototo$pos[length(tototo$hei)])
    if (close < 100) {
      method = "iter"
    }
    else {
      method = "iter2"
    }
  }
  if (is.null(init.thresh)) {
    if (mean(x) > 1000) {
      init.thresh <- quantile(hohoho$hei, 0.95)
    }
    else {
      init.thresh <- median(hohoho$hei)/2.2
    }
  }
  MSE <- function(x, y) {
    X <- cbind(1, x)
    qr.X <- qr(X)
    b <- t(qr.Q(qr.X)) %*% y
    R <- qr.R(qr.X)
    beta <- as.vector(backsolve(R, b))
    fit <- X %*% beta
    res <- list(mse = sum((y - fit)^2), beta = beta[2])
    return(res)
  }
  MSE2 <- function(x, y) {
    X <- cbind(1, x)
    qr.X <- qr(X)
    b <- t(qr.Q(qr.X)) %*% y
    R <- qr.R(qr.X)
    beta <- as.vector(backsolve(R, b))
    fit <- X %*% beta
    mse <- sum((y - fit)^2)
    sst <- sum((y - mean(y))^2)
    r2 <- 1 - (mse/sst)
    res <- list(mse = mse, beta = beta, r2 = r2)
    return(res)
  }
  MSE3 <- function(mix, miy) {
    X <- cbind(1, mix)
    qr.X <- qr(X)
    b <- t(qr.Q(qr.X)) %*% miy
    R <- qr.R(qr.X)
    beta <- as.vector(backsolve(R, b))
    fit <- X %*% beta
    mse <- sum((miy - fit)^2)
    sst <- sum((miy - mean(miy))^2)
    r2 <- 1 - (mse/sst)
    res <- list(mse = mse, beta = beta, r2 = r2)
    return(res)
  }
  thresh = init.thresh
  roxy <- big.peaks.col(x[1:length(x)], thresh)
  if (!is.null(reducing)) {
    nono <- which(roxy$pos %in% reducing)
    roxy <- lapply(roxy, function(x) {
      x[nono]
    })
  }
  nnn <- length(roxy$pos)
  fdsa = 1
  while (nnn < length(ladder)) {
    if (fdsa == 1) {
      cat("\nReducing threshold 2x to find ladder \n")
    }
    thresh = thresh/2
    roxy <- big.peaks.col(x[1:length(x)], thresh)
    nnn <- length(roxy$pos)
    fdsa <- fdsa + 1
  }
  whot <- length(roxy$pos) * 0.2
  what <- which(roxy$hei == max(roxy$hei))
  roxy <- separate(roxy, shift = sep.index, type = "pos")
  if (method == "iter") {
    ii <- which(roxy$hei == max(roxy$hei)) + 1
    iii <- length(roxy$hei)
    roxy <- list(pos = roxy$pos[ii:iii], hei = roxy$hei[ii:iii])
    step1 <- combn(roxy$pos[1:attempt], 3)
    step2 <- apply(step1/10, 2, MSE, y = ladder[1:3])
    mse <- unlist(lapply(step2, function(x) {
      x$mse
    }))
    covs <- apply(step1, 2, function(x, y) {
      cov(x, y)
    }, y = ladder[1:3])
    step2 <- mse * covs
    step3 <- step1[, which(step2 < sort(step2, decreasing = FALSE)[20])]
    step4 <- apply(step3, 2, function(x, y) {
      which(y %in% x)
    }, y = roxy$pos)
    caller <- function(roxy, www, ladder.call, x) {
      threshold <- length(x)
      posi <- numeric()
      fact2 <- length(ladder.call)
      expect <- roxy$pos[www]
      xxx <- ladder.call[c(1:3)]
      modx <- lm(expect ~ poly(xxx, degree = 1))
      expecto <- predict(modx, data.frame(xxx = ladder.call))
      ladder.call <- ladder.call[which(expecto < threshold * 
                                         0.85)]
      available <- length(roxy$pos) - length(ladder.call)
      ava2 <- length(ladder.call) - abs(available)
      if (available > 0) {
        if ((length(ladder.call) - 1) < 3) {
          tope <- length(ladder.call)
        }
        else {
          tope <- length(ladder.call) - 1
        }
      }
      else {
        tope <- ava2 - 2
      }
      expect <- rep(NA, tope + 1)
      for (i in 3:tope) {
        if (i == 3 & i != tope) {
          expect[1:3] <- roxy$pos[www]
          xxx <- ladder.call[c(1:3)]
          mod <- MSE2(xxx, expect[1:3])
          beta <- (mod)$beta
          expecto <- as.vector(beta[1] + matrix(ladder.call) %*% 
                                 beta[-1])
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          expect[i + 1] <- act[good]
          if (mod$r2 < 0.9) {
            i = tope
          }
        }
        if (i > 3 & i <= 5) {
          xx <- ladder.call[c(1:i)]
          mod <- MSE2(xx, expect[1:i])
          beta <- (mod)$beta
          expecto <- as.vector(beta[1] + matrix(ladder.call) %*% 
                                 beta[-1])
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          expect[i + 1] <- act[good]
          if (mod$r2 < 0.9) {
            i = tope
          }
        }
        if (i > 5) {
          xx <- cbind(ladder.call[c(1:i)], ladder.call[c(1:i)]^2, 
                      ladder.call[c(1:i)]^3, ladder.call[c(1:i)]^4)
          mod <- MSE2(xx, expect[1:i])
          beta <- (mod)$beta
          if (length(which(is.na(beta))) > 0) {
            beta[which(is.na(beta))] <- 0
          }
          toto <- cbind(matrix(ladder.call), matrix(ladder.call)^2, 
                        matrix(ladder.call)^3, matrix(ladder.call)^4)
          expecto <- cbind(rep(1, dim(toto)[1]), toto) %*% 
            beta
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          expect[i + 1] <- act[good]
          if (is.na(mod$r2)) {
            mod$r2 <- 0.1
          }
          if (mod$r2 < 0.9) {
            i = tope
          }
        }
        if (i == tope & i != 3) {
          if (i < 5) {
            expect[1:3] <- roxy$pos[www]
            xx <- ladder.call[c(1:i)]
          }
          else {
            xx <- cbind(ladder.call[c(1:i)], ladder.call[c(1:i)]^2, 
                        ladder.call[c(1:i)]^3, ladder.call[c(1:i)]^4)
          }
          mod <- MSE2(xx, expect[1:i])
          beta <- (mod)$beta
          if (length(which(is.na(beta))) > 0) {
            beta[which(is.na(beta))] <- 0
          }
          if (i < 5) {
            toto <- cbind(matrix(ladder.call))
          }
          else {
            toto <- cbind(matrix(ladder.call), matrix(ladder.call)^2, 
                          matrix(ladder.call)^3, matrix(ladder.call)^4)
          }
          expecto <- cbind(rep(1, dim(toto)[1]), toto) %*% 
            beta
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          expect[i + 1] <- act[good]
        }
        if (i == tope & i == 3) {
          expect[1:3] <- roxy$pos[www]
        }
      }
      posi <- expect
      tutu <- abs(length(x) - posi)
      posi <- posi[1:which(tutu == min(tutu, na.rm = TRUE))]
      heii <- roxy$hei[which(roxy$pos %in% posi)]
      fact3 <- length(posi)/fact2
      if (length((posi)) < 6) {
        fact <- summary(lm(ladder.call[1:length(posi)] ~ 
                             poly(posi, degree = length((posi)) - 1)))$r.squared * 
          fact3
      }
      else {
        fact <- summary(lm(ladder.call[1:length(posi)] ~ 
                             poly(posi, degree = 5)))$r.squared * fact3
      }
      roxy2 <- list(pos = posi, hei = heii, wei = ladder.call[1:length(posi)], 
                    corr = abs(cor(ladder.call[1:length(posi)], posi)), 
                    error = fact)
      return(roxy2)
    }
    rt <- apply(data.frame(step4), 2, FUN = caller, roxy = roxy, 
                ladder.call = ladder, x = x)
    corrs3 <- unlist(lapply(rt, function(x) {
      x$error
    }))
    roxy3 <- rt[[which(corrs3 == max(corrs3))]]
    if (draw == TRUE) {
      limi <- sort(roxy3$hei, decreasing = TRUE)
      ################################################################################
      #                                                                              #
      #                                   MODIFIED                                   #
      #                                                                              #
      ################################################################################
      plot_temp <- plot_calibr(x = x, roxy3 = roxy3, limi = limi)
      ################################################################################
      #                                                                              #
      #                                   MODIFIED                                   #
      #                                                                              #
      ################################################################################
    }
    roxy <- roxy3
  }
  if (method == "iter2") {
    last <- length(roxy$pos)
    lld <- length(ladder)
    if ((last - attempt) < 0) {
      roxy$wei <- ladder[1:length(roxy$pos)]
      roxy$corr <- 0
      roxy$error <- 0
      if (draw == TRUE) {
        limi <- sort(roxy$hei, decreasing = TRUE)
        
        ################################################################################
        #                                                                              #
        #                                   MODIFIED                                   #
        #                                                                              #
        ################################################################################
        plot_temp <- plot_calibr(x = x, roxy3 = roxy3, limi = limi)
        ################################################################################
        #                                                                              #
        #                                   MODIFIED                                   #
        #                                                                              #
        ################################################################################
      }
    }
    else {
      step1 <- combn(roxy$pos[last:(last - attempt)], 3)
      step2 <- apply(step1/10, 2, MSE, y = ladder[lld:(lld - 
                                                         2)])
      mse <- unlist(lapply(step2, function(x) {
        x$mse
      }))
      covs <- apply(step1, 2, function(x, y) {
        cov(x, y)
      }, y = ladder[lld:(lld - 2)])
      step2 <- mse
      step3 <- step1[, which(step2 < sort(step2, decreasing = FALSE)[20])]
      step4 <- apply(step3, 2, function(x, y) {
        sort(which(y %in% x), decreasing = TRUE)
      }, y = roxy$pos)
      caller <- function(roxy, www, ladder.call, x) {
        threshold <- length(x)
        last3 <- length(ladder.call):(length(ladder.call) - 
                                        2)
        posi <- numeric()
        fact2 <- length(ladder.call)
        expect <- roxy$pos[www]
        xxx <- ladder.call[last3]
        modx <- lm(expect ~ poly(xxx, degree = 1))
        expecto <- predict(modx, data.frame(xxx = ladder.call))
        available <- length(roxy$pos) - length(ladder.call)
        ava2 <- length(ladder.call) - abs(available)
        if (available < 0) {
          tope <- 3
        }
        else {
          tope <- length(ladder.call)
        }
        expect <- rep(NA, tope)
        lenlad <- length(ladder.call)
        for (i in 3:(tope - 1)) {
          if (i == 3 & i != tope) {
            expect[tope:(tope - 2)] <- roxy$pos[www]
            xxx <- ladder.call[last3]
            mod <- MSE3(mix = xxx, miy = expect[tope:(tope - 
                                                        2)])
            beta <- (mod)$beta
            expecto <- as.vector(beta[1] + matrix(sort(ladder.call, 
                                                       decreasing = TRUE)) %*% beta[-1])
            condo <- sort(expect[which(!is.na(expect))], 
                          decreasing = TRUE)
            act <- sort(roxy$pos[-which(roxy$pos %in% 
                                          condo)], decreasing = TRUE)
            yoyo <- abs(expecto[i + 1] - act)
            good <- which(yoyo == min(yoyo, na.rm = TRUE))
            qwer <- i
            qwer2 <- length(expect) - qwer
            expect[qwer2] <- act[good]
            if (mod$r2 < 0.9) {
              i = tope
            }
          }
          else if (i > 3 & i <= 5) {
            xx <- ladder.call[c(lenlad:(lenlad - (i - 
                                                    1)))]
            mod <- MSE3(mix = xx, miy = expect[tope:qwer2])
            beta <- (mod)$beta
            expecto <- as.vector(beta[1] + matrix(sort(ladder.call, 
                                                       decreasing = TRUE)) %*% beta[-1])
            condo <- sort(expect[which(!is.na(expect))], 
                          decreasing = TRUE)
            act <- sort(roxy$pos[-which(roxy$pos %in% 
                                          condo)], decreasing = TRUE)
            yoyo <- abs(expecto[i + 1] - act)
            good <- which(yoyo == min(yoyo, na.rm = TRUE))
            qwer <- i
            qwer2 <- length(expect) - qwer
            expect[qwer2] <- act[good]
            if (mod$r2 < 0.9) {
              i = tope
            }
          }
          else if (i > 5) {
            ladder.call2 <- sort(ladder.call, decreasing = TRUE)
            expect2 <- sort(expect, decreasing = TRUE)
            xx <- cbind(ladder.call2[c(1:i)])
            mod <- MSE3(mix = xx, miy = expect2)
            beta <- (mod)$beta
            if (length(which(is.na(beta))) > 0) {
              beta[which(is.na(beta))] <- 0
            }
            toto <- cbind(matrix(ladder.call2))
            expecto <- cbind(rep(1, dim(toto)[1]), toto) %*% 
              beta
            condo <- sort(expect[which(!is.na(expect))], 
                          decreasing = TRUE)
            act <- sort(roxy$pos[-which(roxy$pos %in% 
                                          condo)], decreasing = TRUE)
            yoyo <- abs(expecto[i + 1] - act)
            good <- which(yoyo == min(yoyo, na.rm = TRUE))
            qwer <- i
            qwer2 <- length(expect) - qwer
            expect[qwer2] <- act[good]
            if (is.na(mod$r2)) {
              mod$r2 <- 0.1
            }
            if (mod$r2 < 0.9) {
              i = tope
            }
          }
          else if (i == tope & i == 3) {
            expect[1:3] <- roxy$pos[www]
          }
        }
        posi <- expect
        heii <- roxy$hei[which(roxy$pos %in% posi)]
        fact3 <- length(posi)/fact2
        if (length((posi)) < 6) {
          fact <- summary(lm(ladder.call[1:length(posi)] ~ 
                               poly(posi, degree = length((posi)) - 1)))$r.squared * 
            fact3
        }
        else {
          fact <- summary(lm(ladder.call[1:length(posi)] ~ 
                               poly(posi, degree = 5)))$r.squared * fact3
        }
        roxy2 <- list(pos = posi, hei = heii, wei = ladder.call[1:length(posi)], 
                      corr = abs(cor(ladder.call[1:length(posi)], 
                                     posi)), error = fact)
        return(roxy2)
      }
      rt <- apply(data.frame(step4), 2, FUN = caller, roxy = roxy, 
                  ladder.call = ladder, x = x)
      corrs3 <- unlist(lapply(rt, function(x) {
        x$error
      }))
      roxy3 <- rt[[which(corrs3 == max(corrs3))[1]]]
      if (draw == TRUE) {
        limi <- sort(roxy3$hei, decreasing = TRUE)
        
        ################################################################################
        #                                                                              #
        #                                   MODIFIED                                   #
        #                                                                              #
        ################################################################################
        plot_temp <- plot_calibr(x = x, roxy3 = roxy3, limi = limi)
        ################################################################################
        #                                                                              #
        #                                   MODIFIED                                   #
        #                                                                              #
        ################################################################################
      }
      roxy <- roxy3
    }
  }
  return(list(res = roxy, plot = plot_temp))
}





################################################################################
#                                                                              #
#                            MY_LADDER.INFO.ATTACH                             #
#                                                                              #
################################################################################
#' *DESCRIPTION*
#' Fonction modifiée de la fonction Fragman::ladder.info.attach(). Les modifications
#' inclut l'insertion de my_find.ladder() ainsi qu'un bout de code permettant
#' de tracer la courbe de calibration. Tout le reste de la fonction est 
#' identique à l'originale. Les modifications sont entourés de blocs de 
#' commentaires "MODIFIED"

my_ladder.info.attach <- function (stored,
                                   ladder,
                                   channel.ladder = NULL,
                                   method = "iter2",
                                   ladd.init.thresh = NULL,
                                   env = parent.frame(),
                                   draw = TRUE,
                                   attempt = 10)
  
{
  all.names <- names(stored)
  dev = 50
  warn = FALSE
  if (is.null(channel.ladder)) {
    channel.ladder <- dim(stored[[1]])[2]  # nombre de channel
  } else {
    channel.ladder <- channel.ladder
  }
  layout(matrix(1:3, nrow = 3, ncol = 1))
  list.ladders <- lapply(stored, function(x) {
    y <- x[, channel.ladder]  # extraction des informations du channel précisé par .fsa
    return(y)
  })
  for (t in 1:length(list.ladders)) {
    attributes(list.ladders[[t]]) <- list(mycomm = names(list.ladders)[t])
  }
  
  ################################################################################
  #                                                                              #
  #                                   MODIFIED                                   #
  #                                                                              #
  ################################################################################
  # Pics détectés
  res <- lapply(
    list.ladders,  # tous les .fsa de la channel précisée
    function(x) my_find.ladder(x = x, # fonction à appliquer
                               ladder = ladder,
                               dev = dev,
                               warn = warn,
                               method = method,
                               init.thresh = ladd.init.thresh,
                               draw = draw,
                               attempt = attempt)$res
  )
  
  env$list.data.covarrubias <- res  # stockage dans l'environnement de travail
  
  
  
  # Plot detection de pics
  plot_calibration <- lapply(
    list.ladders,  # tous les .fsa de la channel précisée
    function(x) my_find.ladder(x = x, # fonction à appliquer
                               ladder = ladder,
                               dev = dev,
                               warn = warn,
                               method = method,
                               init.thresh = ladd.init.thresh,
                               draw = draw,
                               attempt = attempt)$plot
  )
  
  
  
  
  # Plot la courbe de regression de la calibration
  calib_df <- lapply(res, function(x)
    data.frame(
      Size = x$wei,
      Time = x$pos
    )
  )
  
  # plot_calibration_regression <-
  #   lapply(calib_df, function(x)
  #     ggplot(x, aes(x = Size, y = Time)) +
  #       geom_point(size = 3, color = "black") +
  #       ylim(c(min(x$Time), max(x$Time + x$Time*0.1))) +
  #       geom_line(color = "black") +
  #       geom_text(aes(label = Size), vjust = -2, size = 3.5, color = "darkred", alpha = 1) +
  #       geom_smooth(method = "lm", se = FALSE, color = "darkred", linetype = "dashed", alpha = .8) +
  #       labs(
  #         title = "",
  #         x = "Standard size (bp)",  # LIZ
  #         y = "Migration time"
  #       ) +
  #       theme_classic(base_size = 19)
  #   )
  
  plot_calibration_regression <- 
    lapply(calib_df, function(df) {
      
      tscale <- 1000
      
      # 1. Référence attendue : droite rouge pointillée
      ref_model <- lm(Size ~ Time, data = data.frame(
        Size = sort(df$Size),
        Time = sort(df$Time)
      ))
      ref_df <- data.frame(Time = seq(min(df$Time), max(df$Time), length.out = 200))
      ref_df$Size <- predict(ref_model, newdata = ref_df)
      
      # 2. Courbe réelle verte (régression locale)
      green_model <- loess(Size ~ Time, data = df)
      green_df <- data.frame(Time = seq(min(df$Time), max(df$Time), length.out = 200))
      green_df$Size <- predict(green_model, newdata = green_df)
      
      # 3. Plot
      ggplot() +
        #  Ligne rouge pointillée : référence
        geom_line(data = ref_df, aes(x = Time / tscale, y = Size),
                  color = "#FF4444", linetype = "dashed", linewidth = 0.6) +
        
        #  Courbe verte réelle : suit les points
        geom_line(data = green_df, aes(x = Time / tscale, y = Size),
                  color = "darkgreen", linewidth = 0.7) +
        
        #  Points calibrés
        geom_point(data = df, aes(x = Time / tscale, y = Size),
                   size = 2.5, color = "black") +
        
        #  Étiquettes tailles bp
        geom_text(data = df, aes(x = Time / tscale, y = Size, label = Size),
                  vjust = -1.2, color = "darkgreen", size = 4) +
        
        #  Axes, titre, encadrement
        labs(
          title = "Fonction de calibration - GS500LIZ",
          x = "Temps/1000",
          y = "Taille [Bp]"
        ) +
        theme_classic(base_size = 19) +
        theme(
          panel.grid.major = element_line(color = "grey90", size = 0.3),
          panel.grid.minor = element_line(color = "grey95", size = 0.2),
          axis.line = element_line(color = "black", size = 0.6),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "none"
        )
    })
  
  # Stocker les deux plots (calibration et regression) par fichier .fsa
  list_plot_calibration <- vector(mode = "list", length = length(plot_calibration_regression))
  list_plot_calibration <- lapply(list_plot_calibration, function(x) 
    vector(mode = "list", length = 2))
  
  names(list_plot_calibration) <- names(res)
  
  for(i in 1:length(plot_calibration_regression)){
    list_plot_calibration[[i]][[1]] <- plot_calibration[[i]]
    list_plot_calibration[[i]][[2]] <- plot_calibration_regression[[i]]
  }
  
  
  
  ################################################################################
  #                                                                              #
  #                                   MODIFIED                                   #
  #                                                                              #
  ################################################################################
  
  cat("\nSizing process complete. Information has been stored in the environment for posterior functions.\nFor example to be used by the overview2() or score.markers() functions.")
  correlations <- unlist(lapply(res, function(x) {
    x$corr
  }))
  if (length(which(correlations < 0.92)) > 0) {
    cat(paste("\nWe did not find a good ladder in", length(which(correlations < 
                                                                   0.92)), "sample(s). \nIf you wish to correct it you can try one of the following:\n"))
    cat("\na) The value of ladd.init.thresh might be too low, making noisy peaks too be abundant \n     Solution-- make sure your initial value 'init.thresh' is not below 200 RFUs\nb) You can continue your analysis without worrying for those samples or removing. Identify them as:\n     corro <- unlist(lapply(list.data.covarrubias, function(x){x$corr}))\n     (bad <- which(corro < .9999))\nc) MOST IMPORTANT! you can correct manually the bad samples using the 'ladder.corrector()' function providing the names of the bad samples (below), your ladder, and the information from the 'storing.inds' function, type ?ladder.corrector\n\nNames of the bad sample(s):\n")
  }
  # #  layout(matrix(1, 1, 1))
  # bads <- all.names[which(correlations < 0.92)]
  # if (length(bads) > 0) {
  #   return(bads)
  # }
  # 
  # À la fin de my_ladder.info.attach()
  env$list_plot_calibration <- list_plot_calibration
  # return(list_plot_calibration = list_plot_calibration)
}
