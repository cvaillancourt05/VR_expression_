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

library(data.table)

# -----------
# Paramètres
# -----------

fenetres <- c("10kb", "50kb")
data_dir <- "/lustre09/project/6000443/expression_genes/resultats/variants_rares_associes_aux_outliers"
out_dir <- "/lustre09/project/6000443/expression_genes/resultats/scores_iogc"

# --Versions des scores Z à traiter
versions_z <- c("avec_eQTL_atteints", "avec_eQTL_non_atteints", "sans_eQTL_atteints", "sans_eQTL_non_atteints")

# --Fichier de phénotype (pour distinguer atteints/non-atteints)
fichier_pheno_global <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/GCbroad_plink.txt"

# ----------------------
# Chargement des données
# ----------------------

entete <- colnames(fread("/lustre09/project/6000443/expression_genes/resultats/retrait_effet_eqtl/phenotypes.txt", nrows = 0)) # --Pour avoir FID et IID
df_cohorte <- data.table(FID_IID = entete[-1])
df_cohorte[, c("FID", "IID") := tstrsplit(FID_IID, "_", fixed = TRUE)]
df_cohorte[, FID_IID := NULL]
df_cohorte[, `:=`(FID = as.character(FID), IID = as.character(IID))]

# --Phénotype global : 0=inconnu, 1=non-atteint, 2=atteint
pheno_global <- fread(fichier_pheno_global, header = FALSE, col.names = c("FID", "IID", "Pheno"))
pheno_global[, IID := as.character(IID)]
iid_non_atteints <- pheno_global[Pheno == 1, IID]
  
# ----------------------------------------
# Fonction pour le calcul des scores IOGC
# ----------------------------------------

calculer_score_iogc <- function(fenetre, version_z, cohorte) {

  input_file <- file.path(data_dir, sprintf("variants_candidats_%s_%s.txt", fenetre, version_z))

  # --Lecture des variants candidats
  candidats <- fread(input_file, select = c("variant_id", "probe_id", "symbol", "individu_outlier"))
  setnames(candidats, "individu_outlier", "IID")
  candidats[, IID := as.character(IID)]

  # --Ajout des données non-atteints (ref. n-atteints) pour les atteints (ref. atteints)
  est_ref_atteints <- grepl("atteints$", version_z) & !grepl("non_atteints", version_z)
  if(est_ref_atteints) {
    version_z_non_atteints <- sub("_atteints$", "_non_atteints", version_z)
    input_file_non_atteints <- file.path(data_dir, sprintf("variants_candidats_%s_%s.txt", fenetre, version_z_non_atteints))

    candidats_non_atteints <- fread(input_file_non_atteints, select = c("variant_id", "probe_id", "symbol", "individu_outlier"))
    setnames(candidats_non_atteints, "individu_outlier", "IID")
    candidats_non_atteints[, IID := as.character(IID)]
    candidats_non_atteints <- candidats_non_atteints[IID %in% iid_non_atteints]

    candidats <- rbind(candidats, candidats_non_atteints)
  }

  # --Déduplication
  candidats_unique <- unique(candidats, by = c("IID", "symbol"))

  # --Calcul
  scores_candidats <- candidats_unique[, .(score_iogc = .N), by = .(IID)]

  # --Inclusion de toute la cohorte
  score_final <- merge(scores_candidats, cohorte, by = "IID", all.y = TRUE)
  score_final[is.na(score_iogc), score_iogc := 0]
  setcolorder(score_final, c("FID", "IID", "score_iogc"))

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