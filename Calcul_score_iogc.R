# --------------------------------------------------------------------------------------------
# Calcul_score_iogc.R
# Calcule le burden de variants rares (IOGC)
#
# Entrées :
#  - Fichier des phénotypes pour la strcuture FID/IID
#  - Variants candidats issus de Filtrage_variants_candidats.R 
#
# Sortie :
#  - Scores bruts selon les différentes conditions
# -------------------------------------------------------------------------------------------

library(tidyverse)
library(data.table)

# -----------
# Paramètres
# -----------

fenetres <- c("10kb", "50kb")
data_dir <- "/home/chloev/links/projects/def-bureau/expression_genes/resultats/variants_rares_associes_aux_outliers/"
out_dir <- "/home/chloev/links/projects/def-bureau/chloev/"

#-- Versions des scores Z à traiter
versions_z <- c("avec_eQTL_atteints", "avec_eQTL_non_atteints", "sans_eQTL_atteints", "sans_eQTL_non_atteints")

# ----------------------
# Chargement des données
# ----------------------

df_cohorte <- read_table("/home/chloev/links/projects/def-bureau/chloev/phenotypes.txt", n_max = 0, show_col_types = FALSE) # --Pour avoir FID et IID
df_cohorte <- data.frame(FID_IID = colnames(df_cohorte)[-1], stringsAsFactors = FALSE) %>%
  separate_wider_delim(FID_IID, delim = "_", names = c("FID", "IID")) %>%
  mutate(FID = as.character(FID), IID = as.character(IID))

# ----------------------------------------
# Fonction pour le calcul des scores IOGC
# ----------------------------------------

calculer_score_iogc <- function(fenetre, version_z, cohorte) {

  input_file <- file.path(data_dir, sprintf("variants_candidats_%s_%s.txt", fenetre, version_z))

  # --Lecture des variants candidats
  candidats <- fread(input_file, select = c("variant_id", "probe_id", "symbol", "individu_outlier"))
  setnames(candidats, "individu_outlier", "IID")
  candidats[, IID := as.character(IID)]

  # --Déduplication
  candidats_unique <- unique(candidats, by = c("IID", "symbol"))

  # --Calcul
  scores_candidats <- candidats_unique[, .(score_iogc = .N, liste_variants = paste(unique(variant_id), collapse = ";")), by = .(IID)]

  # --Inclusion de toute la cohorte
  df_score_final <- right_join(scores_candidats, cohorte, by = "IID")
  score_final <- as.data.table(df_score_final)
  score_final[is.na(score_iogc), `:=`(score_iogc = 0, liste_variants = "")]
  setcolorder(score_final, c("FID", "IID", "score_iogc", "liste_variants"))

  return(score_final)

}


# -----------------------------------
# Boucle pour toutes les combinaisons
# -----------------------------------

for (fenetre in fenetres) {
  for(version_z in versions_z) {

    scores_finaux <- calculer_score_iogc(fenetre, version_z, df_cohorte)

    if(!is.null(scores_finaux)) {
      out_file <- file.path(out_dir, sprintf("score_iogc_%s_%s.txt", fenetre, version_z))
      fwrite(scores_finaux, out_file, sep = "\t", quote = FALSE)
    }
  }
}