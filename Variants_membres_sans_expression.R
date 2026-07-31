# -----------------------------------------------------------------------------
# Variants_membres_sans_expression.R
# Identifier et lister les variants chez les membres d'une famille qui ne 
# possèdent pas de données d'expression
#
# Entrées : 
#   - Fichier comprenant les IID, les FID et les phénotypes
#   - Matrice de statut d'expression par famille
#   - Fichier de mapping entre variants et sondes
#
# Sorties : 
#   - Fichier de résultats par fenêtre contenant la liste des variants
# -----------------------------------------------------------------------------

library(data.table)

data_dir      <- "/lustre09/project/6000443/expression_genes/resultats/variants_rares_associes_aux_outliers/porteurs"
candidats_dir <- "/lustre09/project/6000443/expression_genes/resultats/variants_rares_associes_aux_outliers"
pheno_file    <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/GCbroad_plink.txt"
expr_matrix   <- "/lustre09/project/6000443/expression_genes/Adjusted_expression_values.txt"
out_dir       <- "/lustre09/project/6000443/expression_genes/resultats"

fenetres  <- c("10kb", "50kb")
versions_z <- c("avec_eQTL_atteints", "avec_eQTL_non_atteints", "sans_eQTL_atteints", "sans_eQTL_non_atteints")

# -----------------------------------------
# Individus ayant des données d'expression
# -----------------------------------------
entete_expr <- strsplit(readLines(expr_matrix, n = 1), "\t", fixed = TRUE)[[1]]
iids_avec_expr <- as.character(entete_expr[-1])

# ---------------------------------------------------
# Identification des sujets et de leur famille (FID)
# ---------------------------------------------------
pheno <- fread(pheno_file, header = FALSE, col.names = c("FID", "IID", "Diagnostic"))
pheno[, `:=`(FID = as.character(FID), IID = as.character(IID))]
pheno <- unique(pheno, by = "IID")

# -------------------
# Boucle par fenetre 
# -------------------
tous_resultats <- list()

for (fenetre in fenetres) {
  for (version in versions_z) {
  
    fichier_candidats <- file.path(candidats_dir, sprintf("variants_candidats_%s_%s.txt", fenetre, version))
    
    # --Si le fichier de candidats n'existe pas ou est vide
    if (!file.exists(fichier_candidats) || file.size(fichier_candidats) == 0) next
    
    candidats <- fread(fichier_candidats)
    if (nrow(candidats) == 0) next
    
    candidats[, `:=`(variant_id = as.character(variant_id), probe_id = as.character(probe_id),
                     individu_outlier = as.character(individu_outlier))]
    setnames(candidats, "individu_outlier", "IID_outlier") 

    # --Fusion avec typage strict
    candidats <- merge(candidats, pheno[, .(IID, FID)], by.x = "IID_outlier", by.y = "IID", all.x = FALSE)
    setnames(candidats, "FID", "FID_outlier")
    if (nrow(candidats) == 0) next

    # -------------------------------------------------------
    # Filtrage des variants par le variant_id déjà identifié
    # -------------------------------------------------------
    tsv_files <- list.files(data_dir, pattern = sprintf("porteurs_%s_chr.*\\.tsv$", fenetre), full.names = TRUE)
    if (length(tsv_files) == 0) next

    porteurs <- rbindlist(lapply(tsv_files, function(f) {
      if (file.size(f) == 0) return(NULL)
      dt <- fread(f, header = FALSE, sep = "\t", col.names = c("IID", "variant_id", "geno"))
      dt[, IID := sub("^[^_]+_", "", as.character(IID))]
      dt[, variant_id := as.character(variant_id)]
      
      # --Filtrage rapide
      dt[variant_id %chin% candidats$variant_id]
    }), use.names = TRUE, fill = TRUE)
    
    if (is.null(porteurs) || nrow(porteurs) == 0) next
    
    porteurs <- unique(porteurs, by = c("IID", "variant_id"))
    porteurs <- merge(porteurs, pheno[, .(IID, FID)], by = "IID", all.x = FALSE)  

    # -------------------------------------------------------------------
    # Filtrage afin d'obtenir les variants des individus sans expression
    # -------------------------------------------------------------------
    resultat <- merge(candidats, porteurs, by = "variant_id", allow.cartesian = TRUE)
    resultat <- resultat[IID != IID_outlier]
    resultat <- resultat[!(IID %chin% iids_avec_expr)]

    if (nrow(resultat) == 0) next

    export <- unique(resultat[, .(fenetre = fenetre, version_z = version, IID_outlier, variant_id, probe_id,
                                   symbol, Z_outlier, FID, IID_porteur = IID, geno)])
    
    if (nrow(export) > 0) {
      tous_resultats[[paste(fenetre, version)]] = export
    }
  }
}

# -----------
# Sauvegarde
# -----------
if (length(tous_resultats) > 0) {
  export_final <- rbindlist(tous_resultats, use.names = TRUE, fill = TRUE)
} else {
  export_final <- data.table()
}

out_file <- file.path(out_dir, "membres_sans_expr.txt")

if (nrow(export_final) > 0) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(export_final, out_file, sep = "\t", quote = FALSE)
  message(sprintf("-> %s porteurs sans expression au total, sauvegardes dans %s", nrow(export_final), out_file))
} else {
  if (file.exists(out_file)) file.remove(out_file)
  message("Aucun membre sans expression porteur trouve.")
}