################################################################################
#                                                                              #
#        GENOFREE - MODULE DE CALIBRATION AUTOMATIQUE                          #
#                                                                              #
################################################################################
# CONTEXTE :
# Cette fonction permet de charger dynamiquement des fichiers chromatogrammes
# (.fsa), de les calibrer automatiquement à partir d’un fichier de référence (.txt)
# contenant les informations sur les mix, marqueurs et fluorochromes.
# 
# OBJECTIFS PRINCIPAUX :
# - Charger des fichiers .fsa et un fichier de référence .txt
# - Effectuer une calibration initiale automatique avec un seuil par défaut
# - Identifier les calibrations échouées
# - Appliquer une recalibration itérative sur ces fichiers avec différents seuils
# - Visualiser les résultats via tableaux et graphiques interactifs
# - Fournir une interface fluide, thématisée et accessible avec pagination
#
# DÉPENDANCES :
# - Fragman : calibration des fichiers .fsa
# - ggplot2, DT, shinyjs, etc. : visualisation, interactivité, style
# - GENOFREE_fonctions.R : fonctions personnalisées (my_ladder.info.attach)
################################################################################

################################################################################
#                               LIBRAIRIES                                     #
################################################################################
library(shiny)             
library(shinyFiles)        
library(ggplot2)          
library(shinythemes)       
library(tidyverse)         
library(shinyjs)           # Pour afficher/cacher dynamiquement des éléments
library(openxlsx)          
library(rhandsontable)     
library(bslib)            
library(Fragman)           
library(DT)                # Tableaux dynamiques interactifs
library(fs)                # Fonctions liées au système de fichiers

################################################################################
#                               SOURCES                                        #
################################################################################
# Import des fonctions personnalisées depuis un script externe

source("GENOFREE_fonctions.R")

################################################################################
#                               INTERFACE UI                                   #
################################################################################
# Construction de l’interface utilisateur avec onglets (navbarPage)

ui <- navbarPage(
  title = "GENOFREE",
  id = "navbar",
  theme = shinythemes::shinytheme("cerulean"),
  
  header = tagList(
    shinyjs::useShinyjs(),
    tags$head(tags$style(HTML("table.dataTable {
        border-collapse: collapse !important;
        table-layout: auto !important;
      }
      table.dataTable td, table.dataTable th {
        border: 1px solid #ddd !important;
        padding: 5px;
        font-size: 15px !important;
      }
      table.dataTable th {
        font-weight: bold !important;
        white-space: nowrap !important;
      }
      .tab-content .apropos-tab h2 {
        font-size: 28px !important;
      }
      .tab-content .apropos-tab h3 {
        font-size: 24px !important;
      }
      .tab-content .apropos-tab p {
        font-size: 20px !important;
      }
      .shiny-notification {
        position: fixed !important;
        top: 20px !important;
        right: 40px !important;
        width: 500px !important;
        font-size: 16px !important;
        font-weight: bold !important;
        color: #000000 !important;
      }
      #info {
        font-size: 17px;
        font-weight: bold;
        border: 2px solid #0000ff;
        border-radius: 15px;
        background-color: rgba(76, 80, 180, 0.1);
        color: #333333;
        padding: 15px;
        height: 60px;
        width: 100%;
        box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.1);
        overflow: hidden;
        text-align: center;
        line-height: 60px;
        transition: all 0.3s ease;
        margin: auto;
        display: flex;
        align-items: center;
        justify-content: center;
        position: relative;
      }
      #info:hover {
        box-shadow: 0px 8px 16px rgba(0, 0, 0, 0.2);
        background-color: rgba(76, 80, 150, 0.2);
      }
      .navbar-nav li a {
        font-size: 18px !important;
        font-weight: bold !important;
      }
      .navbar-header h1, .navbar-header h2, .navbar-header h3 {
        font-size: 24px !important;
        font-weight: bold !important;
      }")))
  ),
  
  tabPanel(
    title = "Calibration",
    fluidRow(
      column(
        width = 4,
        h4("1. Choisir un dossier contenant des fichiers .fsa"),
        shinyDirButton("fsa_dir", "Parcourir...", "Sélectionner un dossier"),
        br(), br(),
        verbatimTextOutput("selected_fsa_dir"),
        uiOutput("loading_fsa_ui"),
        
        h4("2. Importer le fichier de référence (.txt)"),
        fileInput("ref_txt", NULL, accept = ".txt"),
        br(),
        
        actionButton("launch_calibration", "1. Calibration initiale", class = "btn btn-primary"),
        uiOutput("progress_calibration_ui"),
        br(), br(),
        
        conditionalPanel(
          condition = "output.calibrationDone == true",
          tagList(
            h4("Fichiers détectés dans le dossier sélectionné"),
            tabsetPanel(
              tabPanel("Calibration initiale", DTOutput("results_table_init")),
              tabPanel("Recalibration", DTOutput("results_table_recalib"))
            ),
            br(),
            actionButton("recalibrate_bad", "2. Recalibrer fichiers échoués", class = "btn btn-danger"),
            uiOutput("progress_recalibration_ui")
          )
        )
      ),
      
      column(
        width = 8,
        uiOutput("plot_pagination"),
        conditionalPanel(
          condition = "output.plotVisible == true",
          br(),
          fluidRow(
            column(6, actionButton("prev_page", "⏮️ Précédent")),
            column(6, actionButton("next_page", "⏭️ Suivant"))
          ),
          br(),
          
          selectInput("selected_file", NULL, choices = NULL, label = NULL),
          conditionalPanel(
            condition = "input.selected_file != ''",
            actionButton("show_specific", "Visualiser un fichier spécifique", class = "btn btn-primary"),
            br(), br(),
            plotOutput("calibration_plot", height = "300px"),
            plotOutput("calibration_curve", height = "300px")
          )
        )
      )
    )
  ),
  
  # Définition des onglets - Calibration + à propos 
  tabPanel("Profil génétique", sidebarLayout(
    sidebarPanel(h4("Fonctionnalités à venir...")),
    mainPanel(h4("Données profil génétique à venir..."))
  )),
  
  tabPanel("À PROPOS", div(class = "apropos-tab",
                           h2("À propos de cette application"),
                          ))
)  

################################################################################
#                               SERVER                                         #
################################################################################
# Comporte toute la logique de l'application :
# - Réactions aux événements utilisateur (calibration, sélection de fichier, navigation)
# - Génération des graphiques et tableaux
# - Suivi des états internes : progression, plots, calibration échouée
# - Recalibration automatique avec mise à jour des résultats graphiques et textuels

server <- function(input, output, session) {
  # Définition des racines pour la navigation dans le système de fichiers
  roots <- c(Home = fs::path_home(), "Ordinateur" = "/")
  shinyDirChoose(input, "fsa_dir", roots = roots, session = session)
  
  # Variables réactives pour stocker les états internes
  fsa_path <- reactiveVal(NULL)                       # Chemin du dossier fsa sélectionné
  my_samples <- reactiveVal()                         # Fichiers .fsa importés
  ladder_used <- reactiveVal()                        # Ladder sélectionnée pour calibration
  all_calibrated <- reactiveVal()                     # Données calibrées
  recalibrated_data <- reactiveVal()                  # Données recalibrées
  bad_files_global <- reactiveVal()                   # Fichiers échoués à la calibration
  plot_store <- reactiveVal()                         # Graphiques à afficher
  current_page <- reactiveVal(1)                      # Page courante pour la pagination
  show_specific_plot <- reactiveVal(FALSE)            # Affichage d’un fichier spécifique
  calibration_done <- reactiveVal(FALSE)              # Statut de la calibration
  
  # Indique si la calibration a été effectuée
  output$calibrationDone <- reactive({ calibration_done() })
  outputOptions(output, "calibrationDone", suspendWhenHidden = FALSE)
  
  # Gestion de la sélection du dossier .fsa avec feedback visuel
  observeEvent(input$fsa_dir, {
    output$loading_fsa_ui <- renderUI({ span("Importation en cours...", style = "color: orange;") })
    shinyjs::delay(100, {
      fsa_path(parseDirPath(roots, input$fsa_dir))
      shinyjs::delay(300, {
        output$loading_fsa_ui <- renderUI({ NULL })  # Retrait du message après l'importation
      })
    })
  })
  
  # Affichage du chemin sélectionné
  output$selected_fsa_dir <- renderText({
    req(fsa_path())
    paste("Dossier sélectionné :", fsa_path())
  })
  
  # Stockage des fichiers .fsa dans une variable réactive
  observeEvent(fsa_path(), {
    req(fsa_path())
    my_samples(storing.inds(folder = fsa_path()))
  })
  
  # Lecture du fichier .txt de référence (mix, fluorochromes, marqueurs)
  ref_data <- reactive({
    req(input$ref_txt)
    read.table(input$ref_txt$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  
  # Construction d’une liste de ladder selon mix et fluorochrome
  ladder_list <- reactive({
    df <- ref_data()
    ladder_list <- list()
    combis <- unique(df[, c("fluo", "mix")])
    for (i in seq_len(nrow(combis))) {
      key <- paste0(combis$fluo[i], "_", combis$mix[i])
      ladder_list[[key]] <- sort(unique(df$pb[df$fluo == combis$fluo[i] & df$mix == combis$mix[i]]))
    }
    ladder_list
  })
  
  # Calibration initiale automatique
  observeEvent(input$launch_calibration, {
    output$progress_calibration_ui <- renderUI({ span("Calibration en cours...", style = "color: orange;") })
    shinyjs::delay(100, {
      req(my_samples(), input$ref_txt)
      list.data.covarrubias <<- list()
      ladd <- ladder_list()[["LIZ_mix1"]]
      ladder_used(ladd)
      
      my_ladder.info.attach(
        stored = my_samples(), 
        ladder = ladd, method = "iter2", 
        ladd.init.thresh = NULL, 
        channel.ladder = 5, 
        draw = TRUE
      )
      
      all <- list.data.covarrubias
      
      # Création du tableau des résultats
      result <- data.frame(
        Fichier = names(all),
        Corr = sapply(all, function(x) {
          if (is.null(x$corr) || x$corr < 0.999 || x$corr > 0.99994) {
            return("N/A")
          } else {
            return(round(x$corr, 6))
          }
        }),
        stringsAsFactors = FALSE
      )
      
      result$Statut <- ifelse(result$Corr == "N/A", "❌ Échec", "✅ Calibré")
      result <- result[order(result$Statut, decreasing = TRUE), ]
      
      output$results_table_init <- renderDT({
        datatable(result, options = list(pageLength = 15), rownames = FALSE) %>%
          formatStyle("Statut", target = "row",
                      backgroundColor = styleEqual(c("✅ Calibré", "❌ Échec"), c("#d9fdd3", "#fddcdc")))
      })
      
      updateSelectInput(session, "selected_file", choices = names(all))
      all_calibrated(all)
      bad_files_global(result$Fichier[result$Statut == "❌ Échec"])
      
      # Réorganisation des graphiques : échoués en premier
      plot_order <- c(result$Fichier[result$Statut == "❌ Échec"], result$Fichier[result$Statut == "✅ Calibré"])
      ordered_plots <- list_plot_calibration[plot_order]
      plot_store(ordered_plots)
      
      current_page(1)
      calibration_done(TRUE)
      
      shinyjs::delay(300, {
        output$progress_calibration_ui <- renderUI({ NULL })
      })
    })
  })
  
  # Recalibration automatique des fichiers échoués
  observeEvent(input$recalibrate_bad, {
    output$progress_recalibration_ui <- renderUI({ span("Recalibration en cours...", style = "color: orange;") })
    shinyjs::delay(100, {
      req(my_samples(), all_calibrated())
      
      corr_min <- 0.999
      corr_max <- 0.99994
      thresh_seq <- seq(300, 2000, by = 100)
      ladd <- ladder_used()
      samples <- my_samples()
      
      list.data.recalibrated <<- list()
      
      bad_files <- names(Filter(function(x) {
        is.null(x$corr) || x$corr < corr_min || x$corr > corr_max
      }, all_calibrated()))
      
      for (file in bad_files) {
        success <- FALSE
        i <- 1
        
        while (!success && i <= length(thresh_seq)) {
          result <- my_ladder.info.attach(
            stored = samples[file],
            ladder = ladd,
            method = "iter2",
            ladd.init.thresh = thresh_seq[i],
            channel.ladder = 5,
            draw = TRUE
          )
          
          if (file %in% names(list.data.covarrubias)) {
            corr_value <- list.data.covarrubias[[file]]$corr
            if (!is.null(corr_value) && corr_value >= corr_min && corr_value <= corr_max) {
              success <- TRUE
              list.data.recalibrated[[file]] <- list.data.covarrubias[[file]]
              attr(list.data.recalibrated[[file]], "used_thresh") <- thresh_seq[i]
            } else {
              i <- i + 1
            }
          } else {
            i <- i + 1
          }
        }
      }
      
      # Affichage des résultats recalibrés
      result <- data.frame(
        Fichier = names(list.data.recalibrated),
        Seuil = sapply(names(list.data.recalibrated), function(name) attr(list.data.recalibrated[[name]], "used_thresh")),
        Corr = sapply(list.data.recalibrated, function(x) round(x$corr, 6)),
        stringsAsFactors = FALSE
      )
      result$Statut <- "✅ Recalibré"
      
      output$results_table_recalib <- renderDT({
        datatable(result, options = list(pageLength = 15), rownames = FALSE) %>%
          formatStyle("Statut", target = "row", backgroundColor = styleEqual("✅ Recalibré", "#d9fdd3"))
      })
      
      # Mise à jour des graphiques
      old_plots <- plot_store()
      updated_plots <- old_plots
      for (f in names(list_plot_calibration)) {
        updated_plots[[f]] <- list_plot_calibration[[f]]
      }
      plot_store(updated_plots)
      
      # Mise à jour des choix pour visualisation spécifique
      updated_choices <- unique(c(names(all_calibrated()), names(list.data.recalibrated)))
      updateSelectInput(session, "selected_file", choices = updated_choices)
      
      shinyjs::delay(300, {
        output$progress_recalibration_ui <- renderUI({ NULL })
      })
    })
  })
  
  # Gestion de la pagination des graphiques (2 fichiers par page)
  output$plot_pagination <- renderUI({
    req(plot_store())
    plots <- plot_store()
    files <- names(plots)
    per_page <- 2
    start <- (current_page() - 1) * per_page + 1
    end <- min(start + per_page - 1, length(files))
    
    panels <- lapply(files[start:end], function(f) {
      fluidRow(
        column(6, renderPlot({ plots[[f]][[1]] })),
        column(6, renderPlot({ plots[[f]][[2]] }))
      )
    })
    do.call(tagList, panels)
  })
  
  # Navigation vers page suivante
  observeEvent(input$next_page, {
    total_pages <- ceiling(length(plot_store()) / 2)
    if (current_page() < total_pages) current_page(current_page() + 1)
  })
  
  # Navigation vers page précédente
  observeEvent(input$prev_page, {
    if (current_page() > 1) current_page(current_page() - 1)
  })
  
  # Contrôle de l'affichage conditionnel des graphiques
  output$plotVisible <- reactive({
    !is.null(plot_store()) && length(plot_store()) > 0
  })
  outputOptions(output, "plotVisible", suspendWhenHidden = FALSE)
  
  # Gestion du bouton "Visualiser un fichier spécifique"
  observeEvent(input$show_specific, {
    show_specific_plot(TRUE)
  })
  
  # Affichage du chromatogramme du fichier sélectionné
  output$calibration_plot <- renderPlot({
    req(show_specific_plot(), input$selected_file, plot_store())
    plot_store()[[input$selected_file]][[1]]
  })
  
  # Affichage de la courbe de calibration du fichier sélectionné
  output$calibration_curve <- renderPlot({
    req(show_specific_plot(), input$selected_file, plot_store())
    plot_store()[[input$selected_file]][[2]]
  })
}

################################################################################
#                               LANCEMENT APP                                  #
################################################################################
shinyApp(ui = ui, server = server)
