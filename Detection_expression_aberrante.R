# ------------------------------------------------------------------------------------------------------------
# Detection_expression_aberrante.R
# Calcul des scores Z et identification des outliers d'expression
#
# Entrées : 
#  - Fichier de l'expression ajustée (sans le retrait du eQTL)
#  - Fichier de l'expression ajustée (avec le retrait du eQTL)
#  - Liste des sondes retenues après filtrage
#  - Liste des phénotypes
# Sorties : 
#  - Matrices de scores Z selon les différentes références -- Z_avec/sans_eQTL_ref_(non)_atteints.txt
#  - Listes des individus retirées au nettoyage -- outliers_globaux_avec/sans.txt
# -------------------------------------------------------------------------------------------------------------

library(data.table)
library(matrixStats)
library(readxl)

# ------------------------
# Préparation des données
# ------------------------

# --Phénotype pour identifier atteints et non-atteints
pheno <- read.table("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/GCbroad_plink.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno) <- c("FID", "IID", "Diagnostic")
pheno$IID <- as.character(pheno$IID)
pheno <- pheno[pheno$Diagnostic != 0, ] # --Supprimer les individus au statut inconnu
pheno$Status <- ifelse(pheno$Diagnostic == 2, 1, 0)

# --Filtrage des sondes pour obtenir seulement les sondes exprimées
# sondes_exprimees <- fread("data/liste_sondes_expr_GC.txt", header = FALSE, col.names ="probe_id")$probe_id
sondes_exprimees <- fread("/lustre09/project/6000443/expression_genes/resultats/liste_sondes_exprimees.txt", header = FALSE, col.names = "probe_id")$probe_id

# --CONDITION AVEC EQTL : expression ajustée, lead cis-eQTL non retiré
expr_avec_eqtl <- fread("/lustre09/project/6000443/expression_genes/Adjusted_expression_values.txt", data.table = FALSE)
colnames(expr_avec_eqtl)[1] <- "probe_id"

expr_avec_eqtl <-expr_avec_eqtl[expr_avec_eqtl$probe_id %in% sondes_exprimees, ]

ids_avec_eqtl <- intersect(colnames(expr_avec_eqtl), pheno$IID)
ids_non_atteints_avec <- pheno$IID[pheno$Status == 0 & pheno$IID %in% ids_avec_eqtl]
ids_atteints_avec     <- pheno$IID[pheno$Status == 1 & pheno$IID %in% ids_avec_eqtl]
 
mat_expr_avec <- as.matrix(expr_avec_eqtl[, ids_avec_eqtl])
rownames(mat_expr_avec) <- expr_avec_eqtl$probe_id

# --CONDITION SANS EQTL : expression ajustée, lead cis-eQTL retiré
expr_sans_eqtl <- fread("/lustre09/project/6000443/expression_genes/resultats/Adjusted_expression_values_without_eqtl.txt", data.table = FALSE)
colnames(expr_sans_eqtl)[1] <- "probe_id"

expr_sans_eqtl <- expr_sans_eqtl[expr_sans_eqtl$probe_id %in% sondes_exprimees, ]

ids_sans_eqtl <- intersect(colnames(expr_sans_eqtl), pheno$IID)
ids_non_atteints_sans <- pheno$IID[pheno$Status == 0 & pheno$IID %in% ids_sans_eqtl]
ids_atteints_sans     <- pheno$IID[pheno$Status == 1 & pheno$IID %in% ids_sans_eqtl]

mat_expr_sans <- as.matrix(expr_sans_eqtl[, ids_sans_eqtl])
rownames(mat_expr_sans) <- expr_sans_eqtl$probe_id

# -------------------------------------------------
# Fonction paramétrable pour calculer les scores Z
# -------------------------------------------------

calculer_scores_z <- function(donnees_matrice, ids_reference, ids_cible) {
  
  matrice_param <- donnees_matrice[, ids_reference, drop = FALSE]

  moyenne <- rowMeans(matrice_param, na.rm = TRUE)
  sd <- rowSds(matrice_param, na.rm = TRUE)
  sd[sd == 0 | is.na(sd)] <- NA

  matrice_cible <- donnees_matrice[, ids_cible, drop = FALSE]
  
  z_scores <- sweep(sweep(matrice_cible, 1, moyenne, "-"), 1, sd, "/")

  return(z_scores)
}

# ---------------------------------------------------------------------
# Calcul inital des scores Z (avant le nettoyage des outliers globaux)
# ---------------------------------------------------------------------

# AVEC EQTL
# Référence : non-atteints -> appliqué à TOUS
z_initial_avec_ref_na <- calculer_scores_z(mat_expr_avec, ids_reference = ids_non_atteints_avec, ids_cible = ids_avec_eqtl)

# Référence : atteints -> appliqué aux atteints uniquement
z_initial_avec_ref_a <- calculer_scores_z(mat_expr_avec, ids_reference = ids_atteints_avec, ids_cible = ids_atteints_avec)


# SANS EQTL
# Référence : non-atteints -> appliqué à TOUS
z_initial_sans_ref_na <- calculer_scores_z(mat_expr_sans, ids_reference = ids_non_atteints_sans, ids_cible = ids_sans_eqtl)

# Référence : atteints -> appliqué aux atteints uniquement
z_initial_sans_ref_a <- calculer_scores_z(mat_expr_sans, ids_reference = ids_atteints_sans, ids_cible = ids_atteints_sans)

# -----------------------------------------------------------
# Nettoyage - identification et retrait des outliers globaux
# -----------------------------------------------------------

# Avec eQTL
outliers_par_ind_avec <- colSums(abs(z_initial_avec_ref_na) >= 2, na.rm = TRUE)
seuil_avec <- quantile(outliers_par_ind_avec, 0.75) + 1.5 * IQR(outliers_par_ind_avec)

ids_gardes_avec <- names(outliers_par_ind_avec[outliers_par_ind_avec <= seuil_avec])
ids_gardes_avec <- setdiff(ids_gardes_avec, "541")
ids_atteints_avec_clean <- intersect(ids_atteints_avec, ids_gardes_avec)
ids_non_atteints_avec_clean <- intersect(ids_non_atteints_avec, ids_gardes_avec)

ids_outliers_globaux_avec <- c(names(outliers_par_ind_avec[outliers_par_ind_avec > seuil_avec]), "541")
write.table(ids_outliers_globaux_avec, "/lustre09/project/6000443/expression_genes/resultats/Zscores/outliers_globaux_avec.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Sans eQTL
outliers_par_ind_sans <- colSums(abs(z_initial_sans_ref_na) >= 2, na.rm = TRUE)
seuil_sans <- quantile(outliers_par_ind_sans, 0.75) + 1.5 * IQR(outliers_par_ind_sans)

ids_gardes_sans <- names(outliers_par_ind_sans[outliers_par_ind_sans <= seuil_sans])
ids_gardes_sans <- setdiff(ids_gardes_sans, "541")
ids_atteints_sans_clean <- intersect(ids_atteints_sans, ids_gardes_sans)
ids_non_atteints_sans_clean <- intersect(ids_non_atteints_sans, ids_gardes_sans)

ids_outliers_globaux_sans <- c(names(outliers_par_ind_sans[outliers_par_ind_sans > seuil_sans]), "541")
write.table(ids_outliers_globaux_sans, "/lustre09/project/6000443/expression_genes/resultats/Zscores/outliers_globaux_sans.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ----------------------------------------------------
# Recalcul final des scores Z sur la matrice nettoyée
# ----------------------------------------------------

# AVEC eQTL
# Référence : non-atteints nettoyés -> appliqué à TOUS les individus retenus
z_final_avec_ref_na <- calculer_scores_z(mat_expr_avec, ids_reference = ids_non_atteints_avec_clean, ids_cible = ids_gardes_avec)
write.table(data.frame(probe_id = rownames(z_final_avec_ref_na), z_final_avec_ref_na, check.names = FALSE), "/lustre09/project/6000443/expression_genes/resultats/Zscores/Z_avec_eQTL_ref_non_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Référence : atteints nettoyés -> appliqué aux atteints uniquement
# Ajout des non-atteints au sous-groupe
z_final_avec_ref_a <- calculer_scores_z(mat_expr_avec, ids_reference = ids_atteints_avec_clean, ids_cible = ids_atteints_avec_clean)
z_final_avec_ref_a <- cbind(z_final_avec_ref_a, z_final_avec_ref_na[, ids_non_atteints_avec_clean, drop = FALSE])
write.table(data.frame(probe_id = rownames(z_final_avec_ref_a), z_final_avec_ref_a, check.names = FALSE), "/lustre09/project/6000443/expression_genes/resultats/Zscores/Z_avec_eQTL_ref_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)



# SANS eQTL
# Référence : non-atteints nettoyés -> appliqué à TOUS les individus retenus
z_final_sans_ref_na <- calculer_scores_z(mat_expr_sans, ids_reference = ids_non_atteints_sans_clean, ids_cible = ids_gardes_sans)
write.table(data.frame(probe_id = rownames(z_final_sans_ref_na), z_final_sans_ref_na, check.names = FALSE), "/lustre09/project/6000443/expression_genes/resultats/Zscores/Z_sans_eQTL_ref_non_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Référence : atteints nettoyés -> appliqué aux atteints uniquement
# Ajout des non-atteints au sous-groupe
z_final_sans_ref_a <- calculer_scores_z(mat_expr_sans, ids_reference = ids_atteints_sans_clean, ids_cible = ids_atteints_sans_clean)
z_final_sans_ref_a <- cbind(z_final_sans_ref_a, z_final_sans_ref_na[, ids_non_atteints_sans_clean, drop = FALSE])
write.table(data.frame(probe_id = rownames(z_final_sans_ref_a), z_final_sans_ref_a, check.names = FALSE), "/lustre09/project/6000443/expression_genes/resultats/Zscores/Z_sans_eQTL_ref_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)