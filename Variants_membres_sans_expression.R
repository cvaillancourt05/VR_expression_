# -----------------------------------------------------------------------------
# Variants_membres_sans_expression.R
# Score + liste détaillée des variants pour les membres sans expression
#
# Sorties : 2 fichiers .txt par fenêtre
# -----------------------------------------------------------------------------

library(data.table)

data_dir <- "/home/chloev/links/projects/def-bureau/chloev"
expr_dir <- "/home/chloev/links/projects/def-bureau/expression_genes"

pheno_file <- file.path(data_dir, "GCbroad_plink.txt")
expr_matrix <- file.path(expr_dir, "valeurs_ajustees_lmekin_lame_sans_batch_all_probes_REML_RMA_GCbr_546sujets_transpose.csv")

fam_expr_file <- file.path(expr_dir, "sondes_exprimees_fam_p05_75percent_549sujets.csv")

merge_dir <- file.path(data_dir, "liste_variants/sorties/merge")
out_dir   <- file.path(data_dir, "membres_sans_expr")

# -----------------------
# Chargement des données
# -----------------------

# --Sujets
pheno <- fread(pheno_file, header = FALSE, col.names = c("FID", "IID", "Diagnostic"))
pheno[, `:=`(FID = as.character(FID), IID = as.character(IID))]

expr_cols <- colnames(fread(expr_matrix, sep = ";", nrows = 0)) 
expr_iid  <- expr_cols[-1] 
expr_iid  <- as.character(expr_iid)

# --Statut d'expression par famille et par sonde
fam_expr <- fread(fam_expr_file)
probe_col <- colnames(fam_expr)[1]
fam_ids   <- colnames(fam_expr)[-1]

# Mise en format long (probe_id, FID, expressed)
fam_expr_long <- melt(fam_expr, id.vars = probe_col, 
                      variable.name = "FID", value.name = "expressed")
setnames(fam_expr_long, probe_col, "probe_id")
fam_expr_long[, FID := as.character(FID)]
fam_expr_true <- fam_expr_long[expressed == TRUE]
fam_expr_true[, expressed := NULL] 


# -----------------------------
# Boucle sur les deux fenêtres
# -----------------------------

fenetres <- c("10kb", "50kb")

for (fenetre in fenetres) {
  
  # --Mapping variant -> probe_id
  mapping_file <- file.path(merge_dir, sprintf("variants_dans_regions_%s.bed", fenetre))
  if (!file.exists(mapping_file)) {
    next
  }
  map <- fread(mapping_file, header = FALSE, col.names = c("variant_id", "probe_gene"))
  map[, probe_id := tstrsplit(probe_gene, "__", fixed = TRUE)[[1]]]
  map[, probe_gene := NULL]
  setkey(map, variant_id)
  
  # --Charger tous les porteurs
  porteur_files <- list.files(merge_dir, 
                              pattern = sprintf("porteurs_%s_chr.*\\.tsv$", fenetre),
                              full.names = TRUE)
  if (length(porteur_files) == 0) {
    next
  }
  porteurs <- rbindlist(lapply(porteur_files, fread, 
                               header = FALSE, sep = "\t",
                               col.names = c("IID", "variant_id", "geno")))
  porteurs[, IID := as.character(IID)]
  porteurs[, variant_id := as.character(variant_id)]
  
  # --Joindre avec le mapping pour obtenir probe_id
  porteurs <- merge(porteurs, map, by = "variant_id", all.x = TRUE)
  porteurs <- porteurs[!is.na(probe_id)]
  if (nrow(porteurs) == 0) {
    next
  }
  
  # --Ajouter le FID du porteur
  porteurs <- merge(porteurs, plink[, .(IID, FID)], by = "IID", all.x = TRUE)
  porteurs <- porteurs[!is.na(FID)]
  porteurs[, has_expr := IID %in% expr_iid]
  
  # --Marquer si la famille est exprimée pour cette sonde
  porteurs <- merge(porteurs, fam_expr_true, by = c("probe_id", "FID"), all.x = TRUE)
  porteurs[is.na(expressed), expressed := FALSE]
  
  # --Sélection : pas d'expression mais famille exprimée
  result <- porteurs[has_expr == FALSE & expressed == TRUE,
                     .(variant_id, probe_id, FID, IID, geno)]
  
  # --Supprimer les doublons
  result <- unique(result)
  
  # --Sauvegarde
  if (nrow(result) > 0) {
    out_file <- file.path(out_dir, sprintf("membres_sans_expr_%s.txt", fenetre))
    fwrite(result, out_file, sep = "\t", quote = FALSE)
  } 
}