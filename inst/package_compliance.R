# These first two should really be moved to global.R
# Load dependencies ----
packs <- c("shiny", "tidyverse", "DT", "shinyjs", "openxlsx", "ggrepel", "scales")
packs_false <- packs[-which(packs %in% installed.packages())]
if (length(packs_false) > 0) {
  install.packages(pkgs = packs_false, dependencies = TRUE)
}
if (packageVersion("scales") < 1.2) {
  install.packages("scales")
}
lapply(packs, library, character.only = TRUE)
rm(packs, packs_false)
