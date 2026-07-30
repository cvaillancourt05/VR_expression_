# ------------------------------------------------------------------
# Outliers_porteurs.R
# Pour chaque variant candidat, compte le nombre de porteurs parmi les
# sujets avec expression, afin de vérifier si
# les variants candidats sont systématiquement portés par une seule
# personne
# ------------------------------------------------------------------

library(data.table)

data_dir      <- "/lustre09/project/6000443/expression_genes/resultats/variants_rares_associes_aux_outliers/porteurs"
candidats_dir <- "/lustre09/project/6000443/expression_genes/resultats/variants_rares_associes_aux_outliers"
pheno_file    <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/GCbroad_plink.txt"
expr_matrix   <- "/lustre09/project/6000443/expression_genes/Adjusted_expression_values.txt"
out_dir       <- "/lustre09/project/6000443/expression_genes/resultats"

fenetres  <- c("10kb", "50kb")
version_z <- c("avec_eQTL_atteints", "avec_eQTL_non_atteints", "sans_eQTL_atteints", "sans_eQTL_non_atteints")

# ----------------------------------------------------
# Chargement des données
# ----------------------------------------------------
# --Identification des sujets avec données d'expression
entete_expr <- strsplit(readLines(expr_matrix, n = 1), "\t", fixed = TRUE)[[1]]
iids_avec_expr <- entete_expr[-1]

# Pour avoir les IID
pheno <- fread(pheno_file, header = FALSE, col.names = c("FID", "IID", "Diagnostic"))
pheno[, `:=`(FID = as.character(FID), IID = as.character(IID))]
pheno <- unique(pheno, by = "IID")

# ------------------------------------------------------------------
# Boucle pour les combinaisons fenêtreXversion
# ------------------------------------------------------------------
tous_resultats <- list()

for (fen in fenetres) {
  for (vz in version_z) {

    fichier_candidats <- file.path(candidats_dir, sprintf("variants_candidats_%s_%s.txt", fen, vz))
    if (!file.exists(fichier_candidats)) {
      next
    }

    candidats <- fread(fichier_candidats)
    candidats[, `:=`(variant_id = as.character(variant_id), probe_id = as.character(probe_id),
                      individu_outlier = as.character(individu_outlier))]
    setnames(candidats, "individu_outlier", "IID_outlier")

    if (nrow(candidats) == 0) next

    # -------------------------------------------------------------------
    # Genotypes des porteurs, filtrés par variant_id des candidats et restreints aux sujets avec expression
    # -------------------------------------------------------------------
    tsv_files <- list.files(data_dir, pattern = sprintf("porteurs_%s_chr.*\\.tsv$", fen), full.names = TRUE)
    porteurs <- rbindlist(lapply(tsv_files, function(f) {
      dt <- fread(f, header = FALSE, sep = "\t", col.names = c("IID", "variant_id", "geno"))
      dt[variant_id %chin% candidats$variant_id]
    }))
    porteurs[, IID := sub("^[^_]+_", "", as.character(IID))]  # -- retire le FID colle (format PLINK)
    porteurs <- unique(porteurs, by = c("IID", "variant_id"))
    porteurs <- porteurs[IID %chin% iids_avec_expr]

    # -------------------------------
    # Nombre de porteurs par variant
    # -------------------------------
    n_porteurs <- porteurs[, .(n_porteurs_avec_expression = uniqueN(IID)), by = variant_id]

    export <- merge(candidats[, .(variant_id, probe_id, symbol, IID_outlier, Z_outlier)],
                     n_porteurs, by = "variant_id", all.x = TRUE)
    export[is.na(n_porteurs_avec_expression), n_porteurs_avec_expression := 0]  # -- ne devrait pas arriver (l'outlier est lui-meme porteur), sert de garde-fou
    export[, singleton := n_porteurs_avec_expression == 1]
    export <- unique(export[, .(fenetre = fen, version_z = vz, variant_id, probe_id, symbol,
                                 IID_outlier, n_porteurs_avec_expression, singleton)])
    tous_resultats[[paste(fen, vz)]] <- export
  }
}

# -----------
# Sauvegarde
# -----------
export_final <- rbindlist(tous_resultats, use.names = TRUE, fill = TRUE)

out_file <- file.path(out_dir, "porteurs_par_variant.csv")
fwrite(export_final, out_file, quote = FALSE)
message(sprintf("-> %s%% des variants sont portés uniquement par 1 individu", round(100 * mean(export_final$singleton), 1)))