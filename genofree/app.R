library(shiny)
library(seqinr)
library(ggplot2)
library(shinythemes)
library(tidyverse)
library(shinyjs)
library(openxlsx)
library(rhandsontable)
library(bslib)

# DATA FOR STANDARD
data(gs500liz)


################################################################################
#################  TABLEAU DE REFERENCE DES ALLELES ############################
################################################################################

ref_allele <- read.table("data/REF_MARQUEUR_MARMOTTE.txt", header = T)
ref_allele$se <- as.numeric(ref_allele$se)

# get ci lower and upper to get the range of microsatellite length according to alleles
ref_allele["ci_lower"] <- ref_allele$pb - ref_allele$se
ref_allele["ci_upper"] <- ref_allele$pb + ref_allele$se

# for control (because one unique marker, needs to increase se)
control_pb <- ref_allele[grepl(pattern = "C", ref_allele$marker), "pb"]
ref_allele[grepl(pattern = "C", ref_allele$marker), "ci_lower"] <- control_pb - 1
ref_allele[grepl(pattern = "C", ref_allele$marker), "ci_upper"] <- control_pb + 1

# rename fluo HEX in VIC
ref_allele$fluo[ref_allele$fluo == "HEX"] <- "VIC"
ref_allele$fluo[ref_allele$fluo == "6FAM"] <- "6-FAM"




################################################################################
########################            FONCTIONS           ########################  
################################################################################
# FONCTIONS --------------------------------------------------------
# 1. PLOT CHROMATOGRAMME CALIBRATION
my_plotabif <- function(abifdata = data,
                        mix, # mix loaded
                        n_ticks,  # add control for tick number
                        chanel = 1,
                        tmin = 1/tscale,
                        tmax = abifdata$Data[["SCAN.1"]]/tscale,
                        tscale = 1000,
                        yscale = 1000,
                        xlab = "Taille [Pb]",
                        ylab = " RFU/1000",
                        irange = (tmin * tscale):(tmax * tscale),
                        x = irange/tscale,
                        xlim = NA,
                        chanel.names = c(1:4, 105),
                        DATA = paste("DATA", chanel.names[chanel], sep = "."),
                        y = abifdata$Data[[DATA]][irange]/yscale,
                        ylim = c(0, max(y)),  # change min(y) to 0
                        main = paste0("Fluochrome : ", abifdata$Data[[paste("DyeN", chanel, sep = ".")]]),
                        calibr,
                        ...) {
  suppressMessages({
    
  # Change color according to fluochrome
    # Determine channel-specific color
    channel_colors <- c("#1C86EE", "olivedrab", "darkorange", "#CD5555", "black")
    color <- channel_colors[chanel]
  
  
  # Apply calibration if provided (ensuring calibration only affects x)
  x <- calibr(irange)  # Ensure calibr is applied consistently
  
  # Set xlim if NA
  if (any(is.na(xlim))) {
    xlim <- c(0, max(x))  # Adjust xlim based on calibrated x values
  }
  

  
  # Filtere references according to mix and channel
  if(mix == "mix1"){
    ref_allele_filtered <- ref_allele[ref_allele$mix == "mix1" &
                                        ref_allele$fluo ==  abifdata$Data[[paste("DyeN", chanel, sep = ".")]],]
  } else if(mix == "mix2"){
    ref_allele_filtered <- ref_allele[ref_allele$mix == "mix2" &
                                        ref_allele$fluo ==  abifdata$Data[[paste("DyeN", chanel, sep = ".")]],]
  }
  
  # Filter unique markers
  unique_markers <- ref_allele_filtered %>%
    distinct(marker, .keep_all = TRUE)  # Keep all columns for the first occurrence of each marker
  
  
  # Create data frame for ggplot
  df_gg <- data.frame(x = x, y = y)
  
  # Filter to remove warning messages
  df_gg <- df_gg %>%
    filter(x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2])

  
  # USE STRIPE TRANSPARENCY
  #' (i.e. un allèle sur deux avec une transparence différente)
  n_allele <- length(ref_allele_filtered$allele)
  rep_temp <- rep(x = c(0.15, 0.3), times = n_allele)
  alpha <- rep_temp[1:n_allele]
  
  # PLOT
  plot_gg <- ggplot(df_gg, mapping = aes(x = x, y = y)) +
    # Expected span of length for microsatellite
    geom_rect(
      data = ref_allele_filtered,
      mapping = aes(
        xmin = ci_lower,
        xmax = ci_upper,
        ymin = 0,
        ymax = ylim[2],
        fill = marker,
        group = allele
      ), 
      inherit.aes = FALSE,  # ignorer les esthétiques globales (sinon ne fonctionne pas)
      alpha = alpha
    ) +
    
    # display allele
    geom_text(data = ref_allele_filtered, aes(x = ci_lower,
                                              y = 3, 
                                              label = allele
                                              ), size = 6, angle = 90, hjust = 0.5, vjust = 1) +
    
    # display marker
    geom_label(data = unique_markers, aes(x = ci_lower,
                                              y = ylim[2], 
                                              label = marker,
                                              group = marker
    ), size = 6,
    label.padding = unit(0.25, "lines"),
    fill = "white",
    color = "black"
    ) +
    
    geom_line(color = color) +
    
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
  
  # ADD color for rectangle per marker
  if(chanel == 5){
    color_values <- rep("black", 14)
  } else{
    color_values <- c("#4169E1", "#EE6363", "#7CCD7C",  "#836FFF", "#0F0F0F", "#FFC125")
  }
  
  plot_gg <- plot_gg + scale_fill_manual(
    values = color_values[1:length(unique(ref_allele_filtered$marker))],
    breaks = unique(arrange(ref_allele_filtered, pb)$marker)  # to be in the same order as in datatable() later
  ) 

  })
  
  return(plot_gg)   # Return the plot directly
}



# 2. PLOT CHROMATOGRAMME NON CALIBRE
my_plotabif_non_calibr <- function(abifdata = data,
                                   chanel = 1,
                                   thres,
                                   tmin = 1/tscale,
                                   tmax = abifdata$Data[["SCAN.1"]]/tscale,
                                   tscale = 1000,
                                   yscale = 1000,
                                   xlab = "Temps/1000",
                                   ylab = "RFU/1000",
                                   main = "Chromatogramme non-calibré",
                                   irange = (tmin * tscale):(tmax * tscale),
                                   x = irange/tscale,
                                   xlim = c(tmin, tmax),
                                   chanel.names = c(1:4, 105),
                                   DATA = paste("DATA", chanel.names[chanel], sep = "."),
                                   y = abifdata$Data[[DATA]][irange]/yscale,
                                   ylim = c(0, max(y)),  # change min(y) to 0
                                   ...) {
  suppressMessages({
    # Create data frame for ggplot
    df_gg <- data.frame(x = x, y = y)
    
    # Filter to remove warning messages
    df_gg <- df_gg %>%
      filter(x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2])
    
    
    # APPLICATION OF THE THRESHOLD
    df_gg$y[df_gg$y < thres[1]] <- 0  # lower threshold
    
    
    # PLOT
    plot_gg <- ggplot(df_gg, mapping = aes(x = x, y = y)) +
      geom_line() +
      
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
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        text = element_text(size = 18)
      )
    
  })
  
  return(plot_gg)   # Return the plot directly
}


# 3. PLOT CALIBRATION
plot_calibration <- function(x, y, xseq, calibr, lizok){
  
  # Rescale for plot
  tscale <- 1000
  
  plot_gg <- ggplot() +
    geom_point(mapping = aes(x = x/tscale, y = y)) + 
    geom_line(mapping = aes(x = xseq/tscale, y = calibr(xseq)), color = "darkgreen") +
    geom_text(mapping = aes(x = x / tscale, y = y, label = lizok[!is.na(lizok)]), nudge_x = .2) +
    
    # Use a default ggplot theme
    theme_minimal() +
    
    # Label xlab, ylab and title
    labs(x = "Temps/1000",
         y = "Taille [Bp]",
         title = "Fonction de calibration - GS500LIZ(75-450)") +
    
    # Set y and x limits (ensure scaling works with calibration)
    
    theme(
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      plot.title = element_text(size = 18, margin = margin(b = 6)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.y = element_text(size = 18, margin = margin(0, 9, 0, 0), face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.text.x = element_text(margin = margin(t = 6)),
      axis.text.y = element_text(margin = margin(r = 6)),
      text = element_text(size = 18)
    )
  
  return(plot_gg)
  
}
  


################################################################################
########################       USER INTERFACE           ########################  
################################################################################
# Define the UI of the app
ui <- bslib::page_fillable(
  theme = shinythemes::shinytheme("cerulean"),
  shinyjs::useShinyjs(),  # To hide elements of user interface
  tags$head(
    # Add custom CSS for table cell borders and other customizations
    tags$style(HTML("
      table.dataTable {
        border-collapse: collapse !important;
        table-layout: auto !important;
      }

      table.dataTable td, table.dataTable th {
        border: 1px solid #ddd !important;
        padding: 5px;
        font-size: 15px !important; /* Reduce font size */
      }

      table.dataTable th {
        font-weight: bold !important;
        white-space: nowrap !important;
      }
      
      /* Augmenter la taille du texte pour l'onglet 'À PROPOS' */
      .tab-content .apropos-tab h2 {
        font-size: 28px !important; /* Taille des titres H2 */
      }
      .tab-content .apropos-tab h3 {
        font-size: 24px !important; /* Taille des titres H3 */
      }
      .tab-content .apropos-tab p {
        font-size: 20px !important; /* Taille des paragraphes */
      }

  }
            
      .shiny-notification {
        position: fixed !important;
        top: 20px !important;
        right: 40px !important;
        width: 500px !important;
        font-size: 16px !important;
        font-weight: bold !important;
        color: #00000 !important;
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
        width: 100%; /* Assure that the element spans the column width */
        box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.1);
        overflow: hidden;
        text-align: center; /* Centers text horizontally */
        line-height: 60px; /* Centers text vertically by matching height */
        transition: all 0.3s ease;
        margin: auto; /* Automatically centers horizontally */
        display: flex; /* Flexbox for universal centering */
        align-items: center; /* Centers vertically */
        justify-content: center; /* Centers horizontally */
        position: relative; /* Keeps it within the column */
      }

      #info:hover {
        box-shadow: 0px 8px 16px rgba(0, 0, 0, 0.2);
        background-color: rgba(76, 80, 150, 0.2);
      }
      
      /* Increase font size for navbar tabs */
      .navbar-nav li a {
        font-size: 18px !important; /* Tab text size */
        font-weight: bold !important; /* Make it bold */
      }

      /* Increase font size for navbar panel titles */
      .navbar-header h1, .navbar-header h2, .navbar-header h3 {
        font-size: 24px !important; /* Title text size */
        font-weight: bold !important; /* Make it bold */
      }
    "))
  ),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition = "input.navbar === 'Calibration' || input.navbar === 'Profil génétique'",
      #COMMON WIDGETS BETWEEN NAV BAR
      # File Inputs
      fileInput(inputId = "fsa_file",  
                label = h4("Fichier .fsa"),
                multiple = FALSE,
                accept = ".fsa",
                placeholder = "Sélectionnez le fichier .fsa",
                buttonLabel = "Parcourir"),
      ),
      
      # Conditional UI for "Calibration" tab
      conditionalPanel(
        condition = "input.navbar === 'Calibration'",  # !!!
        
        # Numeric and Slider Inputs
        numericInput(inputId = "tmin",
                     label = h4("Temps minimum"),
                     value = 1.35,
                     min = 1,
                     max = 4,
                     step = 0.05),
        
        sliderInput(inputId = "thres",
                    label = h4("Seuil localisation pics"),
                    value = c(0.4),
                    min = 0.1,
                    max = 2,
                    step = 0.01)
      ),
      
      # Conditional UI for "Profil génétique" tab
      conditionalPanel(
        condition = "input.navbar === 'Profil génétique'",
        
        fileInput(inputId = "excel_file",
                  label = h4("Fichier séquenceur .xlsx"),
                  multiple = FALSE,
                  accept = ".xlsx",
                  placeholder = "Sélectionnez le fichier genoscreen .xlsx",
                  buttonLabel = "Parcourir"),
        
        # Select Input
        selectInput(inputId = "mix_loaded",
                    label = h4("Mix chargé"), 
                    choices = list("Mix 1" = "mix1", "Mix 2" = "mix2"), 
                    selected = 1),
        
        # REFERECENCE ALLELES TAB
        br(),
        h4("Tableau allèles reférence"),
        actionButton(inputId = "ref_allele_button", label = "Afficher / Cacher"),
        
        
        br(),
        br(),
        # Radio Buttons
        radioButtons(inputId = "fluo_choice",
                     label = h4("Fluochrome"),
                     choices = list("FAM" = 1, "VIC" = 2, "NED" = 3, "PET" = 4, "CONTROL" = 5),
                     selected = 1),
        
        # X and Y Limit Sliders
        sliderInput(inputId = "xlim",
                    label = h4("Intervalle de l'axe des x"),
                    value = c(0, 400),
                    min = 0,
                    max = 600,
                    step = 10),
        
        sliderInput(inputId = "ylim",
                    label = h4("Intervalle de l'axe des y"),
                    value = c(0, 15),
                    min = 0,
                    max = 50,
                    step = 1)
      ),
      width = 2
    ),
    
    mainPanel(
      navbarPage("GENOFREE", id = "navbar",
                 tabPanel("À propos",
                          class = "apropos-tab",  # ID pour cibler ce panneau avec CSS
                          h2(strong("GENOFREE")),
                          strong(p("GENOFREE est une application gratuite [en cours de développement] pour la lecture de chromatogrammes et la détermination du profil génétique d'individus basée sur des données microsatellites. Pour toute question : ", span("pierre-alexandre.quittet@cefe.cnrs.fr", style = "color: blue"), ".", style = "font-size: 20px")),
                          p(strong("Lien GitHub du projet :"), "https://github.com/paquittet/GENOFREE/blob/main/genofree"),
                          br(),
                          h2(strong("TUTORIEL ET FONCTIONNALITÉS")),
                          p(h3("Calibration des données")),
                          p(strong("(1)"), "Importer un fichier .fsa 'genoscreen' à l'aide du bouton 'Parcourir' du widget ", span("Fichier .fsa", style = "font-family:monospace")),
                          p(strong("(2)"), "Ajuster les valeurs du widget ", span("Temps minimum", style = "font-family:monospace"),
                            " de manière à ne faire apparaître que les pics d'intérêt (visible dans la figure 'Chromatogramme attendu')"),
                          p(strong("(3)"), "Ajuster les valeurs du widget ", span("Seuil de localisation des pics", style = "font-family:monospace"),
                            " de manière à retirer les pics superflus [attention, cela nécessite parfois un ajustement précis]"),
                          p("La fonction de calibration doit montrer une droite qui fait apparaître les standards : [75, 100, 139, 150, 160, 200, 300, 350, 400, 450]. 
    Attention, parfois, le standard à 75 pb disparaît de la fonction de calibration mais pas du chromatogramme. 
    Pour une calibration correcte, tous les standards cités précédemment doivent apparaître. 
    Dans ce cas, il est nécessaire d'ajuster le ", 
                            span("seuil de localisation des pics", style = "font-family: monospace"),
                            " et de le diminuer légèrement. Cet ajustement doit parfois être très précis (à 0.1 près)."),
                          
                          br(),
                          p(h3("Lecture des profils génétiques")),
                          h4("Lecture des pics et zoom sur le chromatogramme"),
                          p("Une fois le fichier .fsa importé et la calibration correctement réalisée, nous pouvons déterminer le profil génétique. Le chromatogramme calibré est disponible dans l'onglet 'Profil Génétique'. Cette figure
    montre les différents marqueurs (encadrés) et les allèles connus attendus (lignes verticales de couleur). Il est extrêmement important de bien faire correspondre le numéro de mix du widget ", 
                            span("Mix chargé", style = "font-family:monospace"), 
                            " avec celui du fichier .fsa chargé, sinon les positions théoriques des allèles ne seront pas correctes. Il est possible de faire apparaître le tableau de référence des allèles attendus à l'aide du bouton ", 
                            span("Tableau allèle référence", style = "font-family:monospace"), "."),
                          p("A noter également que le choix du fluochrome à l'aide du widget ", span("fluochrome", style = "font-family:monospace"), 
                            " triera automatiquement le tableau de référence ainsi que les marqueurs et allèles attendus sur le chromatogramme."),
                          p("Une fonctionnalité importante de GENOFREE est la possibilité de ", strong("zoomer sur le chromatogramme"), " de manière à lire plus facilement les pics de fluorescence. Pour zoomer, il suffit de maintenir la souris enfoncée sur le graphique
    et de sélectionner la zone d'intérêt. L'application sélectionnera directement la zone y [0; ymax sélectionné] et l'intervalle des x de la sélection, puis fera apparaître le graphique zoomé. Lorsque le graphique zoomé apparaît, il est possible de 
    sélectionner une partie du graphique zoomé. Le petit encadré ", span("Longueur (pb) de la zone sélectionnée", style = "font-family:monospace"), " donnera la valeur en x (c'est-à-dire la longueur en paire de base) correspondant au maximum de fluorescence (c'est-à-dire l'axe y) de la zone sélectionnée, facilitant ainsi 
    la détermination des allèles."),
                          
                          br(),
                          h4("Modification du fichier Excel Genoscreen directement depuis l'application"),
                          p("Enfin, il est possible de charger dans l'application le fichier Excel Genoscreen fourni avec les lectures du séquenceur.", strong("Ce tableau est modifiable directement depuis l'application"), ". A noter que, pour afficher le tableau, il est nécessaire de charger le fichier .fsa correspondant au fichier Genoscreen. En effet, le tableau
    est automatiquement trié par mix et par fichier .fsa pour faciliter la lecture et la modification. Les marqueurs représentés sur les chromatogrammes sont colorés de la même couleur dans le tableau pour faciliter la modification. Il est possible d'enregistrer les modifications.
    A noter que l'enregistrement des modifications inclura les changements dans le fichier Excel complet Genoscreen, et ces modifications spécifiques au fichier .fsa, aux marqueurs, etc., seront affichées en gras et en bleu dans le nouvel Excel Genoscreen. Il est important de bien donner un nom différent de l'original
    au fichier Excel Genoscreen modifié afin de ne pas écraser l'original. Une fois ce fichier Excel 'modifié' créé, il suffit de le charger à chaque lecture de profil génétique (c'est-à-dire à chaque fichier .fsa) pour le compléter petit à petit avec les lectures et vérifications effectuées.")
                 )
                 
                 ,
                 tabPanel("Calibration",
                          fluidRow(
                            uiOutput(outputId = "title_calibration", inline = TRUE),
                            column(5, plotOutput("plot_peak")),
                            column(5, plotOutput("plot_calibration"))
                          ),
                          fluidRow(br()),
                          fluidRow(br()),
                          fluidRow(
                            column(2),
                            column(6, tags$img(src = "standard_expected.png", width = "100%")),
                            column(2)
                          )
                 ),
                 tabPanel("Profil génétique",
                          fluidRow(
                            uiOutput(outputId = "title1", inline = TRUE),
                            column(4, plotOutput("plot_fsa", brush = brushOpts(id = "plot_fsa_brush", resetOnNew = FALSE))),
                            column(4, plotOutput("plot_fsa_zoom", click = "plot_zoom_click", brush = brushOpts(id = "plot_zoom_fsa_brush", resetOnNew = FALSE))),
                            column(1, style = "text-align: center;", uiOutput(outputId = "title_longueur", inline = TRUE),
                                   verbatimTextOutput("info")),
                            column(2, DT::DTOutput(outputId = "table_allele"))
                          ),
                          
                          fluidRow(br()),
                          fluidRow(
                            uiOutput(outputId = "title2_table", inline = TRUE),
                            downloadButton(outputId = "download", label = "Enregistrer sous"),
                            rHandsontableOutput("table")
                          )
                 )
      ),
      width = 10
    )
  )
)













################################################################################
########################            SERVEUR             ########################  
################################################################################
server <- function(input, output) {

  # MAKE REACTIVITY OF TMIN LONGER LATENCY
  tmin_input <- reactive(input$tmin)
  tmin_delayed <- tmin_input |> 
    debounce(300)
  
  
  # A. REACTIVE VALUES  ----------------------------------------------------------
  
  # REACTIVE TITLE 1
  title_calibration <- reactiveVal(NULL)
  
  # Action pour changer l'état
  observeEvent(input$fsa_file, {
    if(!is.null(title_calibration)){
      title_calibration(TRUE)
    }
  })
  
  # Affichage conditionnel du titre en fonction de l'état de `event_triggered`
  output$title_calibration <- renderUI({
    req(input$fsa_file)
    if (title_calibration()) {
      h2(strong("Calibration"))
    }
  })
  
  
  
  # REACTIVE TITLE 1
  title1_reactive <- reactiveVal(NULL)
  
  # Action pour changer l'état
  observeEvent(input$fsa_file, {
    if(!is.null(title1_reactive)){
      title1_reactive(TRUE)
    }
  })
  
  # Affichage conditionnel du titre en fonction de l'état de `event_triggered`
  output$title1 <- renderUI({
    req(input$fsa_file)
    if (title1_reactive()) {
      h2(strong("Lecture du profil génétique"))
    }
  })
  
  
  # REACTIVE TITLE 2 TABLE
  title2_reactive <- reactiveVal(NULL)
  
  # Action pour changer l'état
  observeEvent(list(input$excel_file, input$fsa_file), {
    req(input$excel_file)
    req(input$fsa_file)
    
    if(!is.null(title2_reactive)){
      title2_reactive(TRUE)
    }
  })
  
  # Affichage conditionnel du titre en fonction de l'état de `event_triggered`
  output$title2_table <- renderUI({
    req(input$excel_file)
    req(input$fsa_file)
    
    if (title2_reactive()) {
     h2(strong("Fichier Genoscreen modifiable"))
    }
  })
  
  
  
  ##################################################################
  # Affichage conditionnel du titre "longueur en pb"
  title_longueur <- reactiveVal(NULL)
  
  observeEvent(input$plot_fsa_brush, {
    
    brush <- input$plot_fsa_brush
    
    if(!is.null(brush)){
      title_longueur(TRUE)
    }
  })
  
  output$title_longueur <- renderUI({
    req(input$fsa_file)
    req(input$plot_fsa_brush)
    
    if (title_longueur()) {
      h4("Longueur (Pb) de la zone sélectionnée")
    }
  })
  ##################################################################
  
  
  # REACTIVE FSA
  fsa_data_reactive <- reactiveValues(data = NULL)
  
  # When a file is uploaded, read it and store the data in fsa_data
  observeEvent(input$fsa_file, {
    req(input$fsa_file)
    fsa_data_reactive$data <- read.abif(input$fsa_file$datapath)
    
    
    # Show notification message to the user
    showNotification(ui = h3("Penser à préciser le mix correspondant au fichier .fsa dans <Mix chargé> !"),
                     type = "error", 
                     duration = 10)
  })
  
  
  
  # REACTIVE XLSX
  xlsx_data_reactive <- reactiveValues(data = NULL, workbook = NULL)
  
  # When a file is uploaded, read it and store the data in fsa_data
  observeEvent(list(input$excel_file, input$fsa_file, input$mix_loaded), {
    req(input$excel_file)
    req(input$fsa_file)
    
    
    # Get which sheet of the genoscreen file to read
    if(input$mix_loaded == "mix1"){
      sheet_number = 2
    } else{
      sheet_number = 3
    }
    
    # # get workbook from loaded excel file
    wb <- openxlsx::loadWorkbook(input$excel_file$datapath)
    
    # Read workbook sheet as data.frame
    data <- openxlsx::readWorkbook(xlsxFile = wb,
                                   sheet = sheet_number,  
                                   startRow = 9,  # because data begin at 9th row
                                   colNames = T)
    
    # Filter by .fsa
    data_filtered <- data[data[,1] %in% input$fsa_file$name,]
    
    # Get rows for reconstruction
    index_modification <- rownames(data_filtered)
    
    # RETURN
    xlsx_data_reactive$data_filtered <- data_filtered
    xlsx_data_reactive$data <- data
    xlsx_data_reactive$index_modification <- index_modification
    xlsx_data_reactive$workbook <- wb
    
  })
  
  
  
  # HIDE ACTION BUTTON WHEN NO .XLSX LOADED
    observe({
    if (is.null(input$excel_file)) {
      shinyjs::hide("download")  # Masquer le bouton si aucun fichier n'est chargé
    } else {
      shinyjs::show("download")  # Afficher le bouton si un fichier est chargé
    }
  })
    
    
    # HIDE LONGUEUR PB WHEN NO ZOOM LOADED
    observe({
      brush <- input$plot_fsa_brush
      
      if (is.null(brush)) {
        shinyjs::hide("info")  # Masquer le bouton si aucun fichier n'est chargé
      } else {
        shinyjs::show("info")  # Afficher le bouton si un fichier est chargé
      }
    })
  
    
  # REACTIVE CALIBRAGE FUNCTION
  calibrage_reactive <- reactiveValues(calibrage_function = NULL,
                                    maxis = NULL,
                                    x = NULL,
                                    y = NULL,
                                    xseq = NULL,
                                    lizok = NULL)

  
  observeEvent(list(input$fsa_file, tmin_delayed(), input$thres), {
    req(input$fsa_file)
    req(tmin_delayed())
    req(input$thres)
    
    # Ensure tmin_delayed() has a valid value before proceeding
    if (is.null(tmin_delayed()) ||
        !is.numeric(tmin_delayed()) ||
        tmin_delayed() < 1) {
      # If tmin_delayed() is NULL or NA, return nothing or clear the plot
      return(calibrage_function = NULL,
             maxis = NULL,
             x = NULL,
             y = NULL,
             xseq = NULL,
             lizok = NULL)
    }
    
    # Peak location calculation (channel 5 is the standard)
    maxis <- peakabif(
      fsa_data_reactive$data,
      chanel = 5,
      npeak = 14,  # Standard with 14 peaks
      tmin = tmin_delayed(),
      fig = FALSE,
      thres = input$thres[1]
    )
    
    # If maxis == only NA, alors ne rien renvoyer
    if (all(is.na(maxis$maxis))) {
      return(list(calibrage_function = NULL,
             maxis = NULL,
             x = NULL,
             y = NULL,
             xseq = NULL,
             lizok = NULL))
    }
    
    # Calibration data
    data(gs500liz)
    lizok <- gs500liz$liz  # longueur en Pb des marqueurs
    lizok[!gs500liz$mask1 | !gs500liz$mask2] <- NA  # rpz les fragments incertains qu'on enlève
    lizok <- lizok[-c(1,2)]  # remove extreme
    y <- lizok[!is.na(lizok)]  
    x <- maxis$maxis[!is.na(lizok)]  # localisation temporelle des pics dans notre échantillon standard
    
    
    # Check if x has enough valid points, others wise cancel app
    if (length(x) < 2 || all(is.na(x))) {
      return(list(calibrage_function = NULL, maxis = NULL, x = NULL, y = NULL, xseq = NULL, lizok = NULL))
    }
    
    # Ensure min and max are finite
    xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 255)
    
    # Calibration function
    calibrage_function <- splinefun(x, y)  # local evaluated the expression here in shiny
    
    # Output
    calibrage_reactive$calibrage_function <- calibrage_function
    calibrage_reactive$maxis <- maxis
    calibrage_reactive$x <- x
    calibrage_reactive$y <- y
    calibrage_reactive$xseq <- xseq
    calibrage_reactive$lizok <- lizok
  })
  
  
  
  
  
  
  # B. RENDER ------------------------------------------------------------------
  
  #########################
  # B.1. PLOT CHROMATOGRAME NON CALIBRE
  #########################
  output$plot_peak <- renderPlot({
    
    # Ensure tmin_delayed() has a valid value before proceeding
    if (is.null(tmin_delayed()) || !is.numeric(tmin_delayed()) || tmin_delayed() < 1) {
      # If tmin_delayed() is NULL or NA, return nothing or clear the plot
      return(NULL)
    }
    
    req(fsa_data_reactive$data)  # Ensure the file is loaded
    req(tmin_delayed())
    req(input$thres)
    
    my_plotabif_non_calibr(abifdata = fsa_data_reactive$data,
                           chanel = 5,  # données fluochromes standard
                           tmin = input$tmin,
                           thres = input$thres,
                           ylim = c(0, input$ylim[2])
                           )
    
  })
  
  
  
  
  #########################
  # B.2. PLOT CALIBRATION
  #########################
  output$plot_calibration <- renderPlot({
    
    # Ensure tmin_delayed() has a valid value before proceeding
    if (is.null(tmin_delayed()) || !is.numeric(tmin_delayed()) || tmin_delayed() < 1) {
      # If tmin_delayed() is NULL or NA, return nothing or clear the plot
      return(NULL)
    }
    
    req(fsa_data_reactive$data)  # Ensure the file is loaded
    req(tmin_delayed())
    req(calibrage_reactive$calibrage_function)
    req(calibrage_reactive$x)
    req(calibrage_reactive$y)
    req(calibrage_reactive$xseq)
    
    
    plot_calibration(x = calibrage_reactive$x,
                     y = calibrage_reactive$y,
                     xseq = calibrage_reactive$xseq, 
                     lizok = calibrage_reactive$lizok,
                     calibr = calibrage_reactive$calibrage_function
                     )
  })
  
  
  
  
  #########################
  # B.3. PLOT CHROMATOGRAMME
  #########################
  output$plot_fsa <- renderPlot({
    
    # Ensure tmin_delayed() has a valid value before proceeding
    if (is.null(tmin_delayed()) || !is.numeric(tmin_delayed()) || tmin_delayed() < 1) {
      # If tmin_delayed() is NULL or NA, return nothing or clear the plot
      return(NULL)
    }
    
    
    # Get data from reactive
    req(fsa_data_reactive$data)  # Ensure the file is loaded
    req(tmin_delayed())
    req(calibrage_reactive$calibrage_function)
    
    # PLOT NO ZOOM
    my_plotabif(
      abifdata = fsa_data_reactive$data,
      mix = input$mix_loaded,
      chanel = as.numeric(input$fluo_choice),
      calibr = calibrage_reactive$calibrage_function,
      xlim = c(as.numeric(input$xlim)),
      ylim = c(as.numeric(input$ylim)),
      n_ticks = as.numeric(as.character(input$nticks))
    )  # local evaluated the expression here in shiny
  })
  
  
  
  
  #########################
  # B.4. PLOT ZOOM
  #########################
  # REACTIVE ZOOM RANGE
  ranges <- reactiveValues(x = NULL, y = NULL)  # for zooming
  
  output$plot_fsa_zoom <- renderPlot({
    # Ensure tmin_delayed() has a valid value before proceeding
    if (is.null(tmin_delayed()) || !is.numeric(tmin_delayed()) || tmin_delayed() < 1) {
      # If tmin_delayed() is NULL or NA, return nothing or clear the plot
      return(NULL)
    }
    
    # PRELIMINARY CHECKS
    req(fsa_data_reactive$data)  # Ensure the file is loaded
    req(ranges$x, ranges$y)  # Explicitly react to ranges, fait disparaitre le plot zoomé si aucun brush
    req(tmin_delayed())
    
    # Plotting the chromatogram with zoom applied
    # plot_gg <- 
    my_plotabif(
      abifdata = fsa_data_reactive$data,
      mix = input$mix_loaded,
      chanel = as.numeric(input$fluo_choice),
      calibr = calibrage_reactive$calibrage_function,
      xlim = c(as.numeric(input$xlim)),
      ylim = c(as.numeric(input$ylim)),
      main = "Zone zoomée"
      # xlim = ranges$x,
      # ylim = c(0, ranges$y[2])
    ) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y)  # enable zoom
  })  
  
  
  # ENABLE ZOOM
  observe({
    brush <- input$plot_fsa_brush  
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(0, brush$ymax)  # borne à 0 ymin
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  

  # GET X CORRESPONDING TO MAXIMUM Y OF SELECTED DATA ON ZOOMEDP PLOT
  ranges_max <- reactiveValues(x = NULL, y = NULL, xmax = NULL)  # for zooming
  
  observe({
    # Get data of zoom box on the zoomed plot
    brush_zoom <- input$plot_zoom_fsa_brush  
    
    # Si zoom, alors 
    if (!is.null(brush_zoom)) {
      # Définir les limites de la zone sélectionnée
      ranges_max$x <- c(brush_zoom$xmin, brush_zoom$xmax)
      ranges_max$y <- c(0, brush_zoom$ymax)  # borne ymin fixée à 0
      
      # Extraire les coordonnées du graphique 
      plot_gg_zoom <- my_plotabif(
        abifdata = fsa_data_reactive$data,
        mix = input$mix_loaded,
        chanel = as.numeric(input$fluo_choice),
        calibr = calibrage_reactive$calibrage_function,
        xlim = c(as.numeric(input$xlim)),
        ylim = c(as.numeric(input$ylim)),
        main = "Zone zoomée"
        # xlim = ranges$x,
        # ylim = c(0, ranges$y[2])
      ) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y)  # enable zoom
      
      # Extraire les data ayant permise la construction du graphique
      plot_data <- ggplot_build(plot_gg_zoom)$data[[4]]
      
      # Filtrer les points à l'intérieur de la zone de sélection
      filtered_points <- subset(plot_data, 
                                x >= brush_zoom$xmin & x <= brush_zoom$xmax & 
                                  y >= 0 & y <= brush_zoom$ymax)
      
      # Trouver le x correspondant au maximum de y dans cette zone
      if (nrow(filtered_points) > 0) {
        max_y <- max(filtered_points$y)
        max_x <- filtered_points$x[which.max(filtered_points$y)]
        ranges_max$max_x <- max_x
        
        # Mettre à jour ranges_max avec les valeurs maximales
      } else {
        ranges_max$x <- NULL
        ranges_max$y <- NULL
      }
    } else {
      # Réinitialiser les valeurs lorsque rien n'est sélectionné
      ranges_max$x <- NULL
      ranges_max$y <- NULL
    }
  })
  
  
  

  #########################
  # B.6. COORDINATES OF THE MOUSE
  #########################
  output$info <- renderText({
    req(input$fsa_file)
    req(tmin_delayed())
    
    paste0(
     round(as.numeric(ranges_max$max_x), 2)
    )
  })
  
  
  
  #########################
  # B.5. PLOT TABLE ALLELE
  #########################
  display_ref_allele_button <- reactiveVal(FALSE)
  observeEvent(input$ref_allele_button, {
    display_ref_allele_button(!display_ref_allele_button())
  })
  
  output$table_allele <- DT::renderDT({
    req(fsa_data_reactive$data)  # Ensure the file is loaded
    req(tmin_delayed())
    req(input$fluo_choice)
    req(input$mix_loaded)
    
    # Sort by size
    ref_allele <- 
      ref_allele |> 
      group_by(mix, fluo) |> 
      arrange(pb) |> 
      as.data.frame()
    
    # Sort by mix and fluochrome
    ref_allele <- ref_allele[
      ref_allele$fluo == fsa_data_reactive$data$Data[[paste0("DyeN.", input$fluo_choice)]] &
        ref_allele$mix == input$mix_loaded, c("marker", "allele", "pb")
    ]
    
    # Change name of variables
    colnames(ref_allele) <- c("Marqueur", "Allele", "Longueur (Pb)")
    
    # Get color depending on levels of marqueur
    if(input$fluo_choice == 5){
      color_values <- rep("black", 10)
      
    } else{
      color_values <- c("#4169E1", "#EE6363", "#7CCD7C",  "#836FFF", "#0F0F0F", "#FFC125")[1:length(unique(ref_allele$Marqueur))]
    }
    
    levels_marker <- unique(ref_allele$Marqueur)
    
    
    if(display_ref_allele_button()){
      
      # Render the table with options to disable selection and search
      DT::datatable(
        ref_allele,
        class = list(stripe = F),  # remove row strip
        options = list(
          dom = 't',  # Disable search box and other controls
          paging = F,  # Disable pagination (i.e. table on different part)
          ordering = FALSE,  # Disable column sorting
          columnDefs = list(
            list(targets = 0, visible = FALSE),  # Hide the first column (row numbers)
            list(targets = "_all", className = "dt[-head|-body]-center")  # Center align all cells
          )
        ),
        selection = 'none'  # Disable row selection
      ) |> 
        DT::formatStyle(
          columns = "Marqueur",
          valueColumns = "Marqueur",
          fontWeight = "bold",
          color = DT::styleEqual(levels = levels_marker, values = color_values)
        ) |> 
        DT::formatStyle(
          columns = "Allele",
          valueColumns = "Marqueur",
          fontWeight = "bold",
          color = DT::styleEqual(levels = levels_marker, values = color_values)
        ) |> 
        DT::formatStyle(
          columns = "Longueur (Pb)",
          valueColumns = "Marqueur",
          fontWeight = "bold",
          color = DT::styleEqual(levels = levels_marker, values = color_values)
        )
      
    }
    
  })
  

  ####################
  # EDITABLE DATAFRAME
  ####################
  output$table <- renderRHandsontable({
    req(xlsx_data_reactive$data)
    req(input$excel_file)
    req(input$fsa_file)
    
    # Highlight depending on marker displayed
    ref_fluo <- list("6-FAM" = 1, "VIC" = 2, "NED" = 3, "PET" = 4, "CONTROL" = 5)
    
    # Ne rien faire si fluo == standard
    if(input$fluo_choice != 5){
      # Sort for better readibility
      ref_allele <- 
        ref_allele |> 
        group_by(mix, fluo) |> 
        arrange(pb) |> 
        as.data.frame()
      
      # extract marker corresponding to ongoing fluo
      ref_allele_filtered <- ref_allele[ref_allele$mix %in% input$mix_loaded &
                                          ref_allele$fluo %in% names(ref_fluo[as.numeric(input$fluo_choice)]), ]
      
      ref_allele_filtered <- ref_allele_filtered[!duplicated(ref_allele_filtered$marker),]
      
      marker_visualized <- ref_allele_filtered$marker
      
      color_values <- c("#4169E1", "#EE6363", "#7CCD7C",  "#836FFF", "#0F0F0F", "#FFC125")[1:length(unique(ref_allele_filtered$marker))]
      
      # Convert R objects into JavaScript-friendly JSON
      marker_visualized_json <- jsonlite::toJSON(marker_visualized, auto_unbox = F)
      color_values_json <- jsonlite::toJSON(color_values, auto_unbox = FALSE)
      
      } 
 
    
    # Renders the table and applies JavaScript renderer
    table <- rhandsontable(xlsx_data_reactive$data_filtered, height = 300) |> 
      hot_col(col = 1:ncol(xlsx_data_reactive$data_filtered), valign = "htCenter", halign = "htCenter")
    
    if(input$fluo_choice != 5){
      table <- table |> 
        hot_cols(renderer = paste0(
          "function (instance, td, row, col, prop, value, cellProperties) {",
          "  Handsontable.renderers.TextRenderer.apply(this, arguments);",
          "  if(value && ", marker_visualized_json, ".includes(value)) {",
          "    var color_values = ", color_values_json, ";",  # Inject the color array
          "    var index = ", marker_visualized_json, ".indexOf(value);",  # Get index of the value
          "    var color = color_values[index % color_values.length];",  # Select color based on index
          "    td.style.background = color;",  # Apply the selected color
          "  }",
          "}"
        ))
    }
    
    table
  })
  
  
  
  output$download <- downloadHandler(
    filename = function(){
      paste0("genoscreen_completed_", Sys.Date(), ".xlsx")
    },
    
    content = function(file){
      req(xlsx_data_reactive$workbook)
      req(hot_to_r(input$table)) 
      
      # Define a custom font style
      font_style <- openxlsx::createStyle(
        fontColour = "#0000ff",  # Set font color (blue)
        textDecoration = "bold",  # Bold text
        halign = "center",
        valign = "center",
        border = "TopBottomLeftRight"
      )
      
      # Load workbook
      wb <- xlsx_data_reactive$workbook  # Ensure `xlsx_data_reactive` is properly defined and accessible in your app
      
      # Load rows with modification
      index_modification <- xlsx_data_reactive$index_modification  # Ensure this exists and has valid data
      
      # Load data
      data <- xlsx_data_reactive$data  # Ensure this exists and has valid data
      
      # Add modifications to the data
      modified_data <- hot_to_r(input$table)  # Ensure `hot_to_r(input$table)` is properly returning a valid data frame
      data[index_modification, ] <- modified_data
      
      # Get the correct sheet to write to
      if (input$mix_loaded == "mix1") {
        sheet_number <- 2
      } else {
        sheet_number <- 3
      }
      
      # Write modifications to the workbook
      for (i in index_modification) {
        # Write the data row
        openxlsx::writeData(
          wb = wb,
          sheet = sheet_number,
          x = t(unlist(data[i, ])),  # Write row as a horizontal vector
          startRow = i,
          colNames = FALSE
        )
        
        # Apply the font style to the written data
        openxlsx::addStyle(
          wb = wb,
          sheet = sheet_number,
          style = font_style,
          rows = i,
          cols = 1:ncol(data),  # Apply the style to all columns in the row
          gridExpand = TRUE     # Expand the style across all specified cells
        )
      }
      
      # Save the workbook to the `file` location provided by `downloadHandler`
      openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
      
      # Show notification message to the user
      showNotification(
        ui = h3("Les modifications apportées depuis l'application dans en gras et de couleur bleu dans le nouveau fichier .xlsx"),
        type = "message",
        duration = 5
      )
    }
  )
}  # end serveur




# Run the application 
shinyApp(ui = ui, server = server)
