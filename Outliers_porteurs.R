# ------------------------------------------------------------------
# Outliers_porteurs.R
# Pour chaque variant candidat, compte le nombre de porteurs parmi les
# sujets avec expression, afin de vérifier si
# les variants candidats sont systématiquement portés par une seule
# personne
# ------------------------------------------------------------------

library(data.table)

candidats_dir <- "/lustre09/project/6000443/expression_genes/resultats/variants_rares_associes_aux_outliers"
out_dir       <- "/lustre09/project/6000443/expression_genes/resultats"

fenetres  <- c("10kb", "50kb")
version_z <- c("avec_eQTL_atteints", "avec_eQTL_non_atteints", "sans_eQTL_atteints", "sans_eQTL_non_atteints")

# ------------------------------------------------------------------
# Boucle pour les combinaisons fenêtre x version
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

    if (nrow(candidats) == 0) next

    # -------------------------------------------------------------------
    # Nombre de porteurs (outliers) distincts par variants
    # -------------------------------------------------------------------
    n_porteurs <- candidats[, .(n_porteurs_avec_expression = uniqueN(individu_outlier)), by = variant_id]
 
    export <- merge(candidats[, .(variant_id, probe_id, symbol, individu_outlier, Z_outlier)],
                     n_porteurs, by = "variant_id")
    export[, singleton := n_porteurs_avec_expression == 1]
    export <- unique(export[, .(fenetre = fen, version_z = vz, variant_id, probe_id, symbol,
                                 individu_outlier, n_porteurs_avec_expression, singleton)])
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