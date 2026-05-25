library(data.table)
library(matrixStats)

# ========================
# Préparation des données
# ========================

# Phénotype pour identifier atteints et non-atteints
pheno <- read.table("data/GCbroad_plink.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno) <- c("FID", "IID", "Diagnostic")
pheno$IID <- as.character(pheno$IID)
pheno <- pheno[pheno$Diagnostic != 0, ]
pheno$Status <- ifelse(pheno$Diagnostic == 2, 1, 0)

# Filtrage des sondes pour obtenir seulement les sondes exprimées
#sondes_exprimees <- fread("data/liste_sondes_expr_GC.txt", header = FALSE, col.names ="probe_id")$probe_id
sondes_exprimees <- fread("results/liste_sondes_exprimees.txt", header = FALSE, col.names = "probe_id")$probe_id

# Données ajustées sans retrait du lead cis-eQTL (AVEC)
expr_avec_eqtl <- fread("data/Adjusted_expression_values.txt", data.table = FALSE)
colnames(expr_avec_eqtl)[1] <- "probe_id"

expr_avec_eqtl <-expr_avec_eqtl[expr_avec_eqtl$probe_id %in% sondes_exprimees, ]

ids_avec_eqtl <- intersect(colnames(expr_avec_eqtl), pheno$IID)
ids_non_atteints_avec <- pheno$IID[pheno$Status == 0 & pheno$IID %in% ids_avec_eqtl]
ids_atteints_avec     <- pheno$IID[pheno$Status == 1 & pheno$IID %in% ids_avec_eqtl]
#Matrice 
mat_expr_avec <- as.matrix(expr_avec_eqtl[, ids_avec_eqtl])
rownames(mat_expr_avec) <- expr_avec_eqtl$probe_id

# Données ajustées avec retrait du lead cis-eQTL (SANS)
expr_sans_eqtl <- fread("results/eQTL/Adjusted_expression_values_without_eqtl.txt", data.table = FALSE)
colnames(expr_sans_eqtl)[1:4] <- c("Chr", "start", "end", "probe_id")

expr_sans_eqtl <- expr_sans_eqtl[expr_sans_eqtl$probe_id %in% sondes_exprimees, ]

ids_sans_eqtl <- intersect(colnames(expr_sans_eqtl), pheno$IID)
ids_non_atteints_sans <- pheno$IID[pheno$Status == 0 & pheno$IID %in% ids_sans_eqtl]
ids_atteints_sans     <- pheno$IID[pheno$Status == 1 & pheno$IID %in% ids_sans_eqtl]
# Matrice
mat_expr_sans <- as.matrix(expr_sans_eqtl[, ids_sans_eqtl])
rownames(mat_expr_sans) <- expr_sans_eqtl$probe_id

# =================================================
# Fonction paramétrable pour calculer les scores Z
# =================================================

calculer_scores_z <- function(donnees_matrice, ids_reference, ids_cible) {
  
  matrice_param <- donnees_matrice[, ids_reference, drop = FALSE]

  moyenne <- rowMeans(matrice_param, na.rm = TRUE)
  sd <- rowSds(matrice_param, na.rm = TRUE)
  sd[sd == 0] <- NA

  matrice_cible <- donnees_matrice[, ids_cible, drop = FALSE]
  
  z_scores <- sweep(sweep(matrice_cible, 1, moyenne, "-"), 1, sd, "/")

  return(z_scores)
}

# ===========================
# Calcul inital des scores Z
# ===========================

# AVEC EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_initial_avec_ref_na <- calculer_scores_z(mat_expr_avec, ids_reference = ids_non_atteints_avec, ids_cible = ids_avec_eqtl)

# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_initial_avec_ref_a <- calculer_scores_z(mat_expr_avec, ids_reference = ids_atteints_avec, ids_cible = ids_atteints_avec)


# SANS EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_initial_sans_ref_na <- calculer_scores_z(mat_expr_sans, ids_reference = ids_non_atteints_sans, ids_cible = ids_sans_eqtl)

# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_initial_sans_ref_a <- calculer_scores_z(mat_expr_sans, ids_reference = ids_atteints_sans, ids_cible = ids_atteints_sans)

# ================================================
# Nettoyage - identification des outliers globaux
# ================================================

# Avec eQTL
outliers_par_ind_avec <- colSums(abs(z_initial_avec_ref_na) >= 2, na.rm = TRUE)
ids_gardes_avec <- names(outliers_par_ind_avec[outliers_par_ind_avec <= 100])
ids_outliers_globaux_avec <- names(outliers_par_ind_avec[outliers_par_ind_avec > 100])

ids_atteints_avec_clean <- intersect(ids_atteints_avec, ids_gardes_avec)
ids_non_atteints_avec_clean <- intersect(ids_non_atteints_avec, ids_gardes_avec)

write.table(ids_outliers_globaux_avec, "results/Zscores/outliers_globaux_avec.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Sans eQTL
outliers_par_ind_sans <- colSums(abs(z_initial_sans_ref_na) >= 2, na.rm = TRUE)
ids_gardes_sans <- names(outliers_par_ind_sans[outliers_par_ind_sans <= 100])
ids_outliers_globaux_sans <- names(outliers_par_ind_sans[outliers_par_ind_sans > 100])

ids_atteints_sans_clean <- intersect(ids_atteints_sans, ids_gardes_sans)
ids_non_atteints_sans_clean <- intersect(ids_non_atteints_sans, ids_gardes_sans)

write.table(ids_outliers_globaux_sans, "results/Zscores/outliers_globaux_sans.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# ======================
# Recalcul des scores Z
# ======================

# AVEC EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_final_avec_ref_na <- calculer_scores_z(mat_expr_avec, ids_reference = ids_non_atteints_avec_clean, ids_cible = ids_gardes_avec)
write.table(data.frame(probe_id = rownames(z_final_avec_ref_na), z_final_avec_ref_na, check.names = FALSE), "results/Zscores/Z_avec_eQTL_ref_non_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_final_avec_ref_a <- calculer_scores_z(mat_expr_avec, ids_reference = ids_atteints_avec_clean, ids_cible = ids_atteints_avec_clean)
write.table(data.frame(probe_id = rownames(z_final_avec_ref_a), z_final_avec_ref_a, check.names = FALSE), "results/Zscores/Z_avec_eQTL_ref_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# SANS EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_final_sans_ref_na <- calculer_scores_z(mat_expr_sans, ids_reference = ids_non_atteints_sans_clean, ids_cible = ids_gardes_sans)
write.table(data.frame(probe_id = rownames(z_final_sans_ref_na), z_final_sans_ref_na, check.names = FALSE), "results/Zscores/Z_sans_eQTL_ref_non_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_final_sans_ref_a <- calculer_scores_z(mat_expr_sans, ids_reference = ids_atteints_sans_clean, ids_cible = ids_atteints_sans_clean)
write.table(data.frame(probe_id = rownames(z_final_sans_ref_a), z_final_sans_ref_a, check.names = FALSE), "results/Zscores/Z_sans_eQTL_ref_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)