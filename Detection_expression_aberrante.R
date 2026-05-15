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

# Données ajustées SANS retrait du lead cis-eQTL
expr_avec_eqtl <- fread("data/Adjusted_expression_values.txt", data.table = FALSE)
colnames(expr_avec_eqtl)[1] <- "probe_id"

# Données ajustées AVEC retrait du lead cis-eQTL
expr_sans_eqtl <- fread("results/eQTL/test/Adjusted_expression_values_without_eqtl.txt", data.table = FALSE)
#expr_sans_eqtl <- expr_avec_eqtl #à supprimer plus tard
colnames(expr_sans_eqtl)[1:4] <- c("Chr", "start", "end", "probe_id")

common_ids <- intersect(colnames(expr_avec_eqtl), pheno$IID)
common_ids <- intersect(common_ids, colnames(expr_sans_eqtl))
common_probes <- intersect(expr_avec_eqtl$probe_id, expr_sans_eqtl$probe_id)

ids_non_atteints <- pheno$IID[pheno$Status == 0 & pheno$IID %in% common_ids]
ids_atteints <- pheno$IID[pheno$Status == 1 & pheno$IID %in% common_ids]

# Matrice d'expression avec eQTL
mat_expr_avec <- as.matrix(expr_avec_eqtl[expr_avec_eqtl$probe_id %in% common_probes, common_ids])
rownames(mat_expr_avec) <- common_probes

# Matrice d'expression sans eQTL
mat_expr_sans <- as.matrix(expr_sans_eqtl[expr_sans_eqtl$probe_id %in% common_probes, common_ids])
rownames(mat_expr_sans) <- common_probes

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
z_initial_avec_ref_na <- calculer_scores_z(mat_expr_avec, ids_reference = ids_non_atteints, ids_cible = common_ids)
#write.table(data.frame(probe_id = rownames(z_initial_avec_ref_na), z_initial_avec_ref_na, check.names = FALSE), "results/Zscores/test.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_initial_avec_ref_a <- calculer_scores_z(mat_expr_avec, ids_reference = ids_atteints, ids_cible = ids_atteints)

# SANS EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_initial_sans_ref_na <- calculer_scores_z(mat_expr_sans, ids_reference = ids_non_atteints, ids_cible = common_ids)
# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_initial_sans_ref_a <- calculer_scores_z(mat_expr_sans, ids_reference = ids_atteints, ids_cible = ids_atteints)

# ================================================
# Nettoyage - identification des outliers globaux
# ================================================
#outliers_par_individu <- colSums(abs(z_initial_avec_ref_na) >= 2, na.rm = TRUE)
#ids_gardes <- names(outliers_par_individu[outliers_par_individu <= 100]) #Rejette tous les individus

outliers_par_individu <- colSums(abs(z_initial_avec_ref_na) >= 2, na.rm = TRUE)
seuil_outlier <- mean(outliers_par_individu) + 3 * sd(outliers_par_individu) #Considéré comme un outlier global si son nombre de sondes aberrantes est anormalement élevé par rapport au reste du groupe
ids_gardes <- names(outliers_par_individu[outliers_par_individu <= seuil_outlier])

ids_atteints_clean <- intersect(ids_atteints, ids_gardes)
ids_non_atteints_clean <- intersect(ids_non_atteints, ids_gardes)

# ======================
# Recalcul des scores Z
# ======================

# AVEC EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_final_avec_ref_na <- calculer_scores_z(mat_expr_avec, ids_reference = ids_non_atteints_clean, ids_cible = ids_gardes)
write.table(data.frame(probe_id = rownames(z_final_avec_ref_na), z_final_avec_ref_na, check.names = FALSE), "results/Zscores/Z_avec_eQTL_ref_non_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_final_avec_ref_a <- calculer_scores_z(mat_expr_avec, ids_reference = ids_atteints_clean, ids_cible = ids_atteints_clean)
write.table(data.frame(probe_id = rownames(z_final_avec_ref_a), z_final_avec_ref_a, check.names = FALSE), "results/Zscores/Z_avec_eQTL_ref_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# SANS EQTL
# 1- Centré/Réduit sur les non-atteints initiaux, appliqué à TOUS les individus
z_final_sans_ref_na <- calculer_scores_z(mat_expr_sans, ids_reference = ids_non_atteints_clean, ids_cible = ids_gardes)
write.table(data.frame(probe_id = rownames(z_final_sans_ref_na), z_final_sans_ref_na, check.names = FALSE), "results/Zscores/Z_sans_eQTL_ref_non_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 2- Centré/Réduit sur les atteints initiaux, appliqué aux atteints uniquement
z_final_sans_ref_a <- calculer_scores_z(mat_expr_sans, ids_reference = ids_atteints_clean, ids_cible = ids_atteints_clean)
write.table(data.frame(probe_id = rownames(z_final_sans_ref_a), z_final_sans_ref_a, check.names = FALSE), "results/Zscores/Z_sans_eQTL_ref_atteints.txt", sep = "\t", row.names = FALSE, quote = FALSE)