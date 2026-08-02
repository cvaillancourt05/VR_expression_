# -----------------------------------------------------------------------------------------------------------------------------
# Retrait_effet_eQTL.R
# Retrait de l'effet du lead cis-eQTL sur l'expression de chaque sonde
#
# Entrées : 
#  - Résultats de l'analyse TensorQTL
#  - Fichier de la matrice d'expression ajustée
#  - Fichier au format .raw du génotype
# Sortie : 
#  - Fichier de la matrice d'expression ajustée selon l'effet du lead cis-eQTL -- Adjusted_expression_values_without_eqtl.txt
# -----------------------------------------------------------------------------------------------------------------------------

library(data.table)

# -----------------------
# Chargement des données
# -----------------------

eqtl_results <- fread("/lustre09/project/6000443/expression_genes/resultats/retrait_effet_eqtl/tensorqtl_out.cis_qtl.txt") 
expression_mat <- fread("/lustre09/project/6000443/expression_genes/Adjusted_expression_values.txt")
raw_geno <- fread("/lustre09/project/6000443/expression_genes/resultats/retrait_effet_eqtl/genotypes_eqtl.raw")

setnames(expression_mat, old = 1, new = "probe")

# -------------------------
# Préparation des matrices
# -------------------------

# -- électionner le lead cis-eQTL de chaque sonde (variant avec la plus petite p-value)
setorder(eqtl_results, phenotype_id, pval_nominal)
lead_eqtls <- eqtl_results[!duplicated(phenotype_id)]

# --Identification
raw_geno_samples_clean <- sub("^[^_]+_", "", as.character(raw_geno$IID))
mat_samples <- intersect(colnames(expression_mat)[-1], raw_geno_samples_clean)

expression_clean <- copy(expression_mat)

# --Restreindre la table des lead eQTL aux sondes réellement présentes dans la matrice
lead_eqtls_clean <- lead_eqtls[phenotype_id %in% expression_clean$probe]
setorder(lead_eqtls_clean, phenotype_id)
setorder(expression_clean, probe)

# --Matrice du génotype
colnames(raw_geno) <- gsub("_[A-Za-z0-9]+$", "", colnames(raw_geno))
raw_geno_filtered <- raw_geno[sub("^[^_]+_", "", as.character(IID)) %in% mat_samples]
cols_to_keep <- c("IID", intersect(colnames(raw_geno_filtered), lead_eqtls_clean$variant_id))
geno_clean <- raw_geno_filtered[, ..cols_to_keep]

geno_clean[, IID := sub("^[^_]+_", "", as.character(IID))]
geno_t <- transpose(geno_clean, make.names = "IID", keep.names = "variant_id")

# ----------------------------------------------------
# Retrait de l'effet du lead cis-eQTL sonde par sonde
# ----------------------------------------------------

# --Pour chaque sonde, régresser l'expression sur le génotype du lead eQTL et 
# remplacer les valeurs originales par les résidus recentrés

for (i in 1:nrow(lead_eqtls_clean)) {
  gene <- lead_eqtls_clean$phenotype_id[i]
  variant <- lead_eqtls_clean$variant_id[i]

  if (gene %in% expression_clean$probe && variant %in% geno_t$variant_id) {
    
    y <- as.numeric(expression_clean[probe == gene, ..mat_samples])
    x <- as.numeric(geno_t[variant_id == variant, ..mat_samples])
    
    model <- lm(y ~ x, na.action = na.exclude)
    adjusted_values <- residuals(model) + mean(y, na.rm = TRUE)

    set(expression_clean, 
      i = which(expression_clean$probe == gene), 
      j = mat_samples,
      value = as.list(adjusted_values))
  }
}

# -----------
# Sauvegarde
# -----------

dt_adjusted <- as.data.table(expression_clean, keep.rownames = "probe")
fwrite(expression_clean, "/lustre09/project/6000443/expression_genes/resultats/Adjusted_expression_values_without_eqtl.txt", sep = "\t", quote = FALSE)