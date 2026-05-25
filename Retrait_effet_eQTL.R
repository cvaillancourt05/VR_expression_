library(data.table)

# ======================
# Chargement des données
# ======================

eqtl_results <- fread("/lustre09/project/6000443/expression_genes/resultats/tensorqtl/tensorqtl.cis_qtl.txt") 
expression_mat <- fread("/lustre09/project/6000443/chloev/Adjusted_expression_values.txt")
raw_geno <- fread("/lustre09/project/6000443/chloev/genotypes.raw")

# =========================
# Préparation des matrices
# =========================

lead_eqtls <- eqtl_results[order(pval_nominal), head(.SD, 1), by = phenotype_id]

# Matrice du génotype
rownames_geno <- as.character(raw_geno$IID)
cols_variants <- colnames(raw_geno)[-(1:6)]

mat_geno_t <- as.matrix(raw_geno[, cols_variants, with = FALSE])
rownames(mat_geno_t) <- rownames_geno
mat_geno <- t(mat_geno_t)

# Matrice d'expression
samples <- intersect(colnames(expression_mat)[-1], colnames(mat_geno))
genes_communs <- intersect(expression_mat$probe, lead_eqtls$phenotype_id)
sub_expr <- expression_mat[probe %in% genes_communs]

# Gestion des doublons
if (any(duplicated(sub_expr$probe))) {
  n_dups <- sum(duplicated(sub_expr$probe))
  sub_expr <- sub_expr[!duplicated(probe)]
}

mat_expr <- as.matrix(sub_expr[, ..samples])
rownames(mat_expr) <- sub_expr$probe

# =======================
# Alignement des données
# =======================

lead_eqtls <- lead_eqtls[phenotype_id %in% sub_expr$probe]
lead_eqtls <- lead_eqtls[!duplicated(phenotype_id)]

variants <- lead_eqtls$variant_id
mat_geno <- mat_geno[rownames(mat_geno) %in% variants, samples, drop = FALSE]

# ====================================
# Retrait de l'effet du lead cis-eQTL
# ====================================

# Matrice pour les données ajustées
mat_adjusted_expr <- mat_expr

for (i in 1:nrow(lead_eqtls)) {
  
  gene <- lead_eqtls$phenotype_id[i]
  variant <- lead_eqtls$variant_id[i]

  if (gene %in% rownames(mat_expr) && variant %in% rownames(mat_geno)) {
    
    y <- as.numeric(mat_expr[gene, ])
    x <- as.numeric(mat_geno[variant, ])
    
    model <- lm(y ~ x, na.action = na.exclude)
    mat_adjusted_expr[gene, ] <- residuals(model) + mean(y, na.rm = TRUE)
  }
}

# ===========
# Sauvegarde
# ===========

dt_adjusted <- as.data.table(mat_adjusted_expr, keep.rownames = "probe")
fwrite(dt_adjusted, "/lustre09/project/6000443/expression_genes/resultats/Adjusted_expression_values_without_eqtl.txt", sep = "\t", quote = FALSE)