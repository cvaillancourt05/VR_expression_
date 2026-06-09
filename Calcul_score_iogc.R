# --------------------------------------------------------------------------------------------
# Calcul_score_iogc.R
# Calcule le burden de variants rares (IOGC)
#
# Entrées :
#  - Cohorte pour l'annotation
#  - Variants candidats issus de Filtrage_variants_candidats.R 
#
# Sortie :
#  - Scores bruts selon les différentes conditions --scores_iogc.csv
# -------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(readxl)

# --Chargement des données
cohorte_ids <- read_table("/home/chloev/links/projects/def-bureau/chloev/phenotypes.txt", n_max = 0, show_col_types = FALSE) # --Pour avoir FID et IID
cohorte_ids <- data.frame(FID_IID = colnames(premiere_ligne)[-1], stringsAsFactors = FALSE) %>%
  # Sépare le format FID_IID en deux colonnes distinctes
  separate_wider_delim(FID_IID, delim = "_", names = c("FID", "IID")) %>%
  mutate(FID = as.character(FID), IID = as.character(IID))

prefix = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/variants_rares_associes_aux_outliers/"
fichiers_variants <- c(
  file.path(prefix, "variants_candidats_10kb_avec_eQTL_atteints.txt"),
  file.path(prefix, "variants_candidats_10kb_avec_eQTL_non_atteints.txt"),
  file.path(prefix, "variants_candidats_10kb_sans_eQTL_atteints.txt"),
  file.path(prefix, "variants_candidats_10kb_sans_eQTL_non_atteints.txt"),
  file.path(prefix, "variants_candidats_50kb_avec_eQTL_atteints.txt"),
  file.path(prefix, "variants_candidats_50kb_avec_eQTL_non_atteints.txt"),
  file.path(prefix, "variants_candidats_50kb_sans_eQTL_atteints.txt"),
  file.path(prefix, "variants_candidats_50kb_sans_eQTL_non_atteints.txt")
)

# -----------------------
# Calcul des scores IOGC
# -----------------------

scores_iogc <- cohorte_ids %>% arrange(FID, IID)

for(fichier in fichiers_variants) {

  # --Extraction des noms
  nom_fichier <- basename(fichier)
  nom_condition <- gsub("variants_candidats_|.txt", "", fichier)

  # --Lecture des variants
  variants <- read_delim(fichier, show_col_types = FALSE)

  # --Compte du burden par individu
  iogc_brut <- variants %>% 
    mutate(individu_outlier = as.character(individu_outlier)) %>% 
    distinct(individu_outlier, symbol) %>% 
    count(individu_outlier, name = paste0("IOGC_raw_", nom_condition))

  col_raw <- paste0("IOGC_raw_", nom_condition)
  col_adj <- paste0("IOGC_adj_", nom_condition)

  # --Fusion
  scores_iogc <- scores_iogc %>% 
    left_join(iogc_brut, by = c("IID" = "individu_outlier")) %>% 
    mutate(
      !!col_raw := coalesce(!!sym(col_raw), 0L), 
      !!col_adj := as.vector(scale(!!sym(col_raw)))
    ) 
}

# -----------
# Sauvegarde
# -----------

write_csv(scores_iogc, "/home/chloev/links/projects/def-bureau/chloev/scores_iogc.csv")