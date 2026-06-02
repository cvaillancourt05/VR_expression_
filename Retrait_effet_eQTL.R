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

eqtl_results <- fread("/lustre09/project/6000443/expression_genes/resultats/tensorqtl/tensorqtl.cis_qtl.txt") 
expression_mat <- fread("/lustre09/project/6000443/chloev/Adjusted_expression_values.txt")
raw_geno <- fread("/lustre09/project/6000443/chloev/genotypes.raw")

# -------------------------
# Préparation des matrices
# -------------------------

# --Sélectionner le lead cis-eQTL de chaque sonde
# lead cis-eQTL = variant avec la plus petite p-value
setorder(eqtl_results, phenotype_id, pval_nominal)
lead_eqtls <- eqtl_results[!duplicated(phenotype_id)]

mat_samples <- intersect(colnames(expression_mat)[-1], as.character(raw_geno$IID))
expression_clean <- expression_mat[probe %in% lead_eqtls$phenotype_id] # --resterindre la matrice aux sondes ayant un lead cis-eQTL identifié

# --Garder la sonde associée au variant avec la plus petite p-value
if (any(duplicated(expression_clean$probe))) {
  pvals <- lead_eqtls[, .(phenotype_id, pval_nominal)]
  expression_clean <- merge(expression_clean, pvals, by.x = "probe", by.y = "phenotype_id", all.x = TRUE)
  
  setorder(expression_clean, probe, pval_nominal)
  expression_clean <- expression_clean[!duplicated(probe)]
  expression_clean[, pval_nominal := NULL]
}

# --Aligner les deux tables
lead_eqtls_clean <- lead_eqtls[phenotype_id %in% expression_clean$probe]
setorder(lead_eqtls_clean, phenotype_id)
setorder(expression_clean, probe)


# --Construire la matrice de génotypes transposée
colnames(raw_geno) <- gsub("_[A-Z0-9]$", "", colnames(raw_geno))
variants_presents <- intersect(colnames(raw_geno), lead_eqtls_clean$variant_id)
raw_geno_filtered <- raw_geno[IID %in% mat_samples]
cols_to_keep <- c("IID", intersect(colnames(raw_geno), lead_eqtls$variant_id))
geno_clean <- raw_geno_filtered[, ..cols_to_keep]
geno_t <- transpose(geno_clean, make.names = "IID", keep.names = "variant_id")

# ----------------------------------------------------
# Retrait de l'effet du lead cis-eQTL sonde par sonde
# ----------------------------------------------------

# --Pour chaque sonde, régresser l'expression sur le génotype du lead eQTL et 
# remplacer les valeurs originales par les résidus recentrés

for (i in 1:nrow(lead_eqtls_clean)) {
  gene <- lead_eqtls$phenotype_id[i]
  variant <- lead_eqtls$variant_id[i]

  if (gene %in% rownames(expression_clean) && variant %in% geno_t$variant_id) {
    
    y <- as.numeric(expression_clean[probe == gene, ..mat_samples])
    x <- as.numeric(geno_t[variant_id == variant, ..mat_samples])
    
    model <- lm(y ~ x, na.action = na.exclude)
    ajusted_values <- residuals(model) + mean(y, na.rm = TRUE)

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
fwrite(dt_adjusted, "/lustre09/project/6000443/expression_genes/resultats/Adjusted_expression_values_without_eqtl.txt", sep = "\t", quote = FALSE)