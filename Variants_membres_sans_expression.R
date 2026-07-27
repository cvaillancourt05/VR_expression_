# -----------------------------------------------------------------------------
# Variants_membres_sans_expression.R
# Identifier et lister les variants chez les membres d'une famille qui ne 
# possèdent pas de données d'expression
#
# Entrées : 
#   - Fichier comprenant les IID, les FID et les phénotypes
#   - Matrice de statut d'expression par famille
#   - Fichier de mapping entre variants et sondes
#   - Liste des variants
#
# Sorties : 
#   - Fichier de résultats par fenêtre contenant la liste des variants
# -----------------------------------------------------------------------------

library(data.table)

data_dir <- "/home/chloev/links/projects/def-bureau/chloev"
expr_dir <- "/home/chloev/links/projects/def-bureau/expression_genes"

pheno_file    <- file.path(data_dir, "GCbroad_plink.txt")
expr_matrix   <- file.path(data_dir, "Adjusted_expression_values.txt")
fam_expr_file <- file.path(expr_dir, "sondes_exprimees_fam_p05_75percent_549sujets.csv")
liste_file    <- file.path(data_dir, "liste_variants.csv")
merge_dir     <- file.path(data_dir, "liste_variants/sorties/merge")
out_dir       <- file.path(data_dir, "membres_sans_expr")

# -----------------------------------------
# IID mesurés dans la matrice d'expression
# -----------------------------------------
entete_expr <- strsplit(readLines(expr_matrix, n = 1), "\t", fixed = TRUE)[[1]]
iids_avec_expr <- entete_expr[-1]

# --------------------------
# Identification des sujets
# --------------------------
pheno <- fread(pheno_file, header = FALSE, col.names = c("FID", "IID", "Diagnostic"))
pheno[, `:=`(FID = as.character(FID), IID = as.character(IID))]
pheno <- unique(pheno, by = "IID")

fam_expr <- fread(fam_expr_file)
probe_col <- colnames(fam_expr)[1]
fam_expr_long <- melt(fam_expr, id.vars = probe_col, variable.name = "FID", value.name = "expressed")
setnames(fam_expr_long, probe_col, "probe_id")
fam_expr_long[, FID := as.character(FID)]
fam_expr_true <- fam_expr_long[expressed == TRUE, .(probe_id, FID)]
fam_expr_true[, probe_id := sub("^X(?=[0-9])", "", probe_id, perl = TRUE)]  # --retire le prefixe X artefact
fam_expr_true[, probe_id := as.character(probe_id)]
fam_expr_true <- unique(fam_expr_true, by = c("probe_id", "FID"))
rm(fam_expr, fam_expr_long); gc()

# --------------------------------------------------------------------
# Variants candidats inclus dans les scores IOGC (liste_variants.csv)
# --------------------------------------------------------------------
liste_variants <- fread(liste_file, colClasses = c(IID = "character"))
setnames(liste_variants, "IID", "IID_outlier")

candidats_long <- merge(liste_variants, pheno[, .(IID, FID)], by.x = "IID_outlier", by.y = "IID", all.x = FALSE)
setnames(candidats_long, "FID", "FID_outlier")

# -------------------
# Boucle par fenetre 
# -------------------
fenetres <- c("10kb", "50kb")

for (fen in fenetres) {
  message(sprintf("Fenetre : %s", fen))

  cand_f <- candidats_long[fenetre == fen]
  if (nrow(cand_f) == 0) next

  porteur_files <- list.files(merge_dir, pattern = sprintf("porteurs_%s_chr.*\\.tsv$", fen), full.names = TRUE)
  if (length(porteur_files) == 0) next

  porteurs <- rbindlist(lapply(porteur_files, function(f) {
    dt <- fread(f, header = FALSE, sep = "\t", col.names = c("IID", "variant_id", "geno"))
    dt[variant_id %chin% cand_f$variant_id]           
  }))
  porteurs[, IID := sub("^[0-9]+_", "", as.character(IID))]  
  porteurs <- unique(porteurs, by = c("IID", "variant_id"))
  porteurs <- merge(porteurs, pheno[, .(IID, FID)], by = "IID", all.x = FALSE)

  # --porteurs du meme variant, meme probe_id exact, meme famille
  resultat <- merge(cand_f, porteurs, by = "variant_id", allow.cartesian = TRUE)
  resultat <- resultat[FID == FID_outlier & IID != IID_outlier]

  # --porteurs sans expression
  sans_expr <- resultat[!(IID %chin% iids_avec_expr)]
  sans_expr <- merge(sans_expr, fam_expr_true, by = c("probe_id", "FID"), all.x = FALSE)
  export_sans_expr <- unique(sans_expr[, .(fenetre, version_z, IID_outlier, variant_id, probe_id, FID,
                                            IID_porteur = IID, geno)])

  if (nrow(export_sans_expr) > 0) {
    fwrite(export_sans_expr, file.path(out_dir, sprintf("membres_sans_expr_%s.txt", fen)), sep = "\t", quote = FALSE)
  } else {
    message("Aucun membre sans expression porteur trouvé pour cette fenetre.")
  }
}


