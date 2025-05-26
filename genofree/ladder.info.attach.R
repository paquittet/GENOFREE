library(Fragman)

# 1. Define the folder containing .fsa files 
fsa_folder <- "~/Desktop/Genofree/genofree/fichier .fsa test/Re Fichier fsa test/fsa"

# 2. Read all .fsa files in the folder 
my.samples <- storing.inds(folder = fsa_folder)

# 3. Define the standard ladder sizes (GS500 LIZ)
LIZ <- c(75, 100, 139, 150, 160, 200, 300, 350, 400, 450)

calib_env <- new.env()


# 4. ladder.info.attach
peakdetect <- ladder.info.attach(
  stored = my.samples,
  ladder = LIZ,
  method = "iter2",
  ladd.init.thresh = 1500,
  channel.ladder = 5,
  draw = TRUE,
  prog = TRUE,
  env = calib_env
)