ref_file <- "/Users/medardabiona/Desktop/Genofree/data/REF_MARQUEUR_MARMOTTE.txt"

ref_marmottes <- read.table(ref_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Créer une structure prête à l’emploi : ladder_list[["NED_mix1"]] = c(positions)
ladder_list <- list()

# Identifier toutes les combinaisons fluo + mix
combis <- unique(ref_marmottes[, c("fluo", "mix")])


for (i in seq_len(nrow(combis))) {
  fluo <- combis$fluo[i]
  mix <- combis$mix[i]
  key <- paste0(fluo, "_", mix)
  
  subset <- ref_marmottes[ref_marmottes$fluo == fluo & ref_marmottes$mix == mix, ]
  ladder_list[[key]] <- unique(sort(subset$pb))
}


fluo_input <- "LIZ"
mix_input <- "mix1"
key <- paste0(fluo_input, "_", mix_input)
ladder_used <- ladder_list[[key]]

to.correct <- my_ladder.info.attach(
  stored = my_samples,
  ladder = ladder_used,
  method = "iter2",
  ladd.init.thresh = NULL,
  channel.ladder = 5,
  draw = TRUE
)


