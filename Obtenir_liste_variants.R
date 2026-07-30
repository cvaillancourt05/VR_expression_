# -----------------------------------------------------------
# Obtenir_liste_variants.R
# Lister les variants contribuant au score IOGC par individu
# -----------------------------------------------------------

library(data.table)

dir_candidats <- "/home/chloev/links/projects/def-bureau/expression_genes/resultats/variants_rares_associes_aux_outliers"

fenetres <- c("10kb", "50kb")
versions_z <- c("avec_eQTL_atteints", "avec_eQTL_non_atteints", 
                "sans_eQTL_atteints", "sans_eQTL_non_atteints")

liste_brute <- list()

for (fenetre in fenetres) {
  for (version in versions_z) {
    fichier <- file.path(dir_candidats, sprintf("variants_candidats_%s_%s.txt", fenetre, version))
    if (!file.exists(fichier)) next

    # --Chargement des données utiles
    dt <- fread(fichier, select = c("individu_outlier", "variant_id", "symbol"))
    setnames(dt, "individu_outlier", "IID")
    
    # --Application du filtre du score
    dt_unique_gene <- unique(dt, by = c("IID", "symbol"))
    
    # --Extraction de la colonne variant_id
    liste_brute[[paste(fenetre, version, sep = "_")]] <- dt_unique_gene[, .(variant_id)]
  }
}

# --Fusion de toutes les conditions et retrait des variants en double 
final_global <- rbindlist(liste_brute)
variants_uniques_iogc <- unique(final_global, by = "variant_id")

# --Sauvegarde
fwrite(variants_uniques_iogc, 
       "/lustre09/project/6000443/expression_genes/resultats/liste_variants.txt", 
       col.names = TRUE, 
       sep = "\t")