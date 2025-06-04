# ref_file <- "/Users/medardabiona/Desktop/Genofree/data/REF_MARQUEUR_MARMOTTE.txt"
ref_file <- "/data/REF_MARQUEUR_MARMOTTE.txt"

ref_marmottes <- read.table(ref_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Créer une structure prête à l’emploi : ladder_list[["NED_mix1"]] = c(positions)
ladder_list <- list()

# Identifier toutes les combinaisons fluo + mix
combis <- unique(ref_marmottes[, c("fluo", "mix")])

################################################################################
#                                                                              #
# Ce que fait cette boucle :                                                   #
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




