################################################################################
#                                                                              #
#         PRÉPARATION DES ÉCHELLES DE CALIBRATION PAR FLUO ET MIX              #
#                                                                              #
################################################################################
#' *CONTEXTE GÉNÉRAL* :
#' Ce bloc de code prépare automatiquement toutes les échelles de calibration 
#' (vecteurs de positions en paires de bases) associées aux combinaisons uniques 
#' de fluorochrome (`fluo`) et de mix (`mix`) présentes dans un fichier .txt 
#' de référence. Ces échelles sont ensuite directement exploitables pour calibrer 
#' les chromatogrammes via la fonction `my_ladder.info.attach()`.
#'
#' *OBJECTIFS* :
#' - Lire le fichier `.txt` de référence des marqueurs.
#' - Extraire pour chaque couple `fluochrome + mix` la liste des tailles attendues (`pb`).
#' - Construire un objet `ladder_list` qui associe chaque combinaison à son vecteur d’échelle.
#' - Permettre à l'utilisateur de sélectionner dynamiquement un `fluochrome` et un `mix`.
#'
#' *ENTRÉES* :
#' - `ref_file` : chemin vers le fichier .txt contenant les colonnes `fluo`, `mix`, `pb`, etc.
#'
#' *SORTIES* :
#' - `ladder_list` : une liste nommée, où chaque entrée est du type `"NED_mix1"` ou `"PET_mix2"`
#'                   et contient le vecteur trié et unique des positions attendues (`pb`).
#'

ref_file <- "/Users/medardabiona/Desktop/Genofree/data/REF_MARQUEUR_MARMOTTE.txt"
# ref_file <- "/data/REF_MARQUEUR_MARMOTTE.txt"

ref_marmottes <- read.table(ref_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Créer une structure prête à l’emploi : ladder_list[["NED_mix1"]] = c(positions)
ladder_list <- list()

# Identifier toutes les combinaisons fluo + mix
combis <- unique(ref_marmottes[, c("fluo", "mix")])


################################################################################
#                                                                              #
# Ce que fait la boucle  for:                                                  #
#                                                                              #
# 1. Parcourt chaque combinaison fluo + mix                                    #       
# 2. Construit une clé d’identification sous la forme "FAM_mix1", "LIZ_mix2".  #
# 3. Filtre ref_marmottes pour ne garder que les lignes correspondant à cette  #
#    combinaison                                                               #
# 4. Extrait la colonne pb (= tailles des fragments attendus)                  #
# 5. Trie et déduplique les valeurs → pour créer la ladder attendue pour cette #
#    combinaison                                                               #
# 6. Stocke le vecteur de positions dans ladder_list[["key"]]                  #
#                                                                              #
################################################################################

for (i in seq_len(nrow(combis))) {
  fluo <- combis$fluo[i]
  mix <- combis$mix[i]
  key <- paste0(fluo, "_", mix)
  
  subset <- ref_marmottes[ref_marmottes$fluo == fluo & ref_marmottes$mix == mix, ]
  ladder_list[[key]] <- unique(sort(subset$pb))
}

# Exemple input 
fluo_input <- "LIZ"
mix_input <- "mix1"
key <- paste0(fluo_input, "_", mix_input)
ladder_used <- ladder_list[[key]]




