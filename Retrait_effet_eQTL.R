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

eqtl_results <- fread("/lustre09/project/6000443/expression_genes/resultats/tensorqtl/tensorqtl_out.cis_qtl.txt") 
expression_mat <- fread("/lustre09/project/6000443/chloev/Adjusted_expression_values.txt")
raw_geno <- fread("/home/chloev/links/projects/def-bureau/chloev/retrait_effet_eqtl/genotypes_eqtl.raw")

setnames(expression_mat, old = 1, new = "probe")

# -------------------------
# Préparation des matrices
# -------------------------

# --Sélectionner le lead cis-eQTL de chaque sonde
# lead cis-eQTL = variant avec la plus petite p-value
setorder(eqtl_results, phenotype_id, pval_nominal)
lead_eqtls <- eqtl_results[!duplicated(phenotype_id)]

mat_samples <- intersect(colnames(expression_mat)[-1], sub("^[^_]+_", "", as.character(raw_geno$IID)))
expression_clean <- copy(expression_mat) # --conserver toutes les sondes; seules celles avec un lead eQTL seront ajustées

# --Restreindre la table des lead eQTL aux sondes réellement présentes dans la matrice
lead_eqtls_clean <- lead_eqtls[phenotype_id %in% expression_clean$probe]
setorder(lead_eqtls_clean, phenotype_id)
setorder(expression_clean, probe)


# --Construire la matrice de génotypes transposée
raw_geno[, IID := sub("^[^_]+_", "", as.character(IID))]
colnames(raw_geno) <- gsub("_[ACGT]+$", "", colnames(raw_geno))

# --Fonction pour normaliser un variant_id en clé insensible à l'ordre REF/ALT
make_key <- function(id) {
  parts <- tstrsplit(id, ":", fixed = TRUE)
  chr <- parts[[1]]; pos <- parts[[2]]
  a1 <- parts[[3]]; a2 <- parts[[4]]
  paste(chr, pos, pmin(a1, a2), pmax(a1, a2), sep = ":")
}

# --Construire les clés des deux côtés
lead_eqtls_clean[, match_key := make_key(variant_id)]

raw_ids <- colnames(raw_geno)[colnames(raw_geno) %like% "^chr"]  # exclure FID/IID/PAT/etc.
key_to_rawid <- setNames(raw_ids, make_key(raw_ids))

# --Vérifier le gain avant de continuer
sum(lead_eqtls_clean$match_key %in% names(key_to_rawid))

raw_geno_filtered <- raw_geno[IID %in% mat_samples]
cols_to_keep <- c("IID", intersect(colnames(raw_geno), lead_eqtls_clean$variant_id))
geno_clean <- raw_geno_filtered[, ..cols_to_keep]
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