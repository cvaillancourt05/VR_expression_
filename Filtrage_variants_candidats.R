# --------------------------------------------------------------------------------------------
# Filtrage_variants_candidats.R
# Identifie les variants rares candidats associés aux outliers d'expression
#
# Entrées :
#  - Paires variant-sonde -- variants_dans_regions_<fenetre>.bed
#  - Fréquences alléliques (pour annotation) -- freq_variants_<fenetre>.frq
#  - Matrice de scores Z -- Z_<version>.txt
#
# Sortie :
#  - Liste des variants candidats par fenetre -- variants_candidats_<fenetre>_<version_z>.txt
# --------------------------------------------------------------------------------------------

library(data.table)

# -----------
# Paramètres
# -----------

fenetres <- c("10kb", "50kb")
seuil_z <- 2

data_dir <- "/home/chloev/links/projects/def-bureau/chloev/liste_variants/sorties/merge"
filtre_dir <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
scores_z <- c("avec_eQTL_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_avec_eQTL_ref_atteints.txt", 
            "avec_eQTL_non_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_avec_eQTL_ref_non_atteints.txt",
            "sans_eQTL_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_sans_eQTL_ref_atteints.txt",
            "sans_eQTL_non_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_sans_eQTL_ref_non_atteints.txt")

# -----------------------------------
# Fonction pour filtrer les variants
# -----------------------------------

filtrer_variants <- function(fenetre, version_z, zscore_file) {

  # -----
  # Associations variant-sonde
  # -----
  regions_dt <- fread(file.path(data_dir, sprintf("variants_dans_regions_%s.bed", fenetre)), header = FALSE, col.names = c("variant_id", "probe_gene"))
  regions_dt[, c("probe_id", "symbol") := tstrsplit(probe_gene, "__", fixed = TRUE)]
  regions_dt[, probe_gene := NULL]

  variants_intersect <- unique(regions_dt$variant_id) # --Liste de tous les variants dans les régions cis

  # -----
  # Chargement des génotypes des porteurs
  # -----
  tsv_files <- list.files(data_dir, pattern = sprintf("porteurs_%s_chr.*\\.tsv$", fenetre), full.names = TRUE)
  
  # --Concaténation de tous les chromosomes en une seule table
  geno_porteurs <- rbindlist(lapply(tsv_files, fread, header = FALSE, sep = "\t", col.names = c("IID", "variant_id", "genotype")))
  # --Extraction de l'IID
  geno_porteurs[, IID := sub(".*_", "", IID)]
  geno_porteurs[, IID := as.character(IID)]
  geno_porteurs[, variant_id := as.character(variant_id)]
  
  # --Exclusion de l'individu aberrant
  geno_porteurs <- geno_porteurs[IID != "541"]

  # --On garde que les porteurs de variants dans une région cis
  ids_valides <- regions_dt[["variant_id"]]
  geno_porteurs <- geno_porteurs[variant_id %in% ids_valides]
  gc()

  # -----
  # Jointure variant-sondes x porteurs
  # -----
  setkey(regions_dt, variant_id)
  setkey(geno_porteurs, variant_id)
  work_initial <- merge(regions_dt, geno_porteurs, by = "variant_id", allow.cartesian = TRUE)
  rm(regions_dt, geno_porteurs)
  gc()

  # --Liste des sondes et individus concernés
  genes_cibles <- unique(work_initial[["probe_id"]])
  individus_cibles <- unique(work_initial[["IID"]])

  # -----
  # Chargement et filtrage des scores Z
  # -----
  zscores_entete <- colnames(fread(zscore_file, nrows = 1))
  cols_z <- c(zscores_entete[1], intersect(zscores_entete, individus_cibles))
  zscores <- fread(zscore_file, select = cols_z)
  setnames(zscores, 1, "probe_id")

  sondes_avec_z <- unique(zscores[["probe_id"]])
  work_sans_z <- work_initial[!probe_id %in% sondes_avec_z]
  work_initial <- work_initial[probe_id %in% sondes_avec_z]

  zscores <- zscores[probe_id %in% unique(work_initial[["probe_id"]])]
  z_long <- melt(zscores, id.vars = "probe_id", variable.name = "IID", value.name = "Z")

  rm(zscores)
  gc()

  # --Calcul du statut outlier et de la direction de déviation
  z_long[, IID := as.character(IID)]
  z_long[, is_outlier := abs(Z) >= seuil_z]
  z_long[, direction := fifelse(Z >= seuil_z, "UP", fifelse(Z <= -seuil_z, "DOWN", "none"))]

  # -----
  # Jointure porteurs x scores Z
  # -----
  work_initial[, probe_id := as.character(probe_id)]
  work_initial[, IID := as.character(IID)]

  z_long[, probe_id := as.character(probe_id)]
  setkey(work_initial, probe_id, IID)
  z_sub <- z_long[, .(probe_id, IID, Z, is_outlier, direction)]
  setkey(z_sub, probe_id, IID)

  work <- z_sub[work_initial, on = .(probe_id, IID)]
  # --Les porteurs sans Z sont marqués non-outlier
  work[is.na(is_outlier), `:=`(is_outlier = FALSE, direction = "none", Z = NA)]

  rm(work_initial, z_long, z_sub)
  gc()

  # -----
  # Critère de présence exclusive chez les outliers
  # -----
  # Variant porté par >= 1 outlier et par aucun non-outlier pour la même sonde
  resume <- work[!is.na(is_outlier), .(n_outliers = sum(is_outlier), 
    n_non_outliers = sum(!is_outlier), directions = paste(unique(direction[is_outlier]), 
    collapse = ",")), by = .(variant_id, probe_id, symbol)]
  
  critere_pres <- resume[n_outliers >= 1 & n_non_outliers == 0]
  rm(resume)

  # -----
  # Critère de direction cohérente de la déviation
  # -----
  # Tous les porteurs outliers dévient dans le même sens
  dir_coherente <- work[is_outlier == TRUE, 
    .(dir_sonde = if(all(direction == "UP")) "UP" 
    else if (all(direction == "DOWN")) "DOWN" else "MIXED"), 
    by = .(variant_id, probe_id)]

  critere_dir <- merge(critere_pres, dir_coherente, by = c("variant_id", "probe_id"))
  critere_dir <- critere_dir[dir_sonde %in% c("UP", "DOWN")] # --On exclut MIXED
  rm(critere_pres, dir_coherente)

  # -----
  # Collecte des ID et scores Z des individus outliers porteurs 
  # -----
  work[, combo_key := paste(variant_id, probe_id)]
  critere_dir[, combo_key := paste(variant_id, probe_id)]

  cles_valides <- critere_dir$combo_key
  work_sub <- work[is_outlier == TRUE & combo_key %in% cles_valides]
  rm(work)
  gc()

  outliers_porteurs <- work_sub[, .(variant_id, probe_id, IID, Z = round(Z, 3))]

  # -----
  # Annotation avec les fréquences alléliques
  # -----
  # Les fichiers .frq contiennent déjà uniquement des variants rares; aucun filtre n'est appliqué
  freq <- fread(file.path(data_dir, sprintf("freq_variants_%s.frq", fenetre)))
  setnames(freq, old = c("SNP", "MAF_A", "MAF_U"), 
            new = c("variant_id", "MAF_Atteints", "MAF_Non_atteints"), skip_absent= TRUE)

  dt_freq <- freq[, .(variant_id, MAF_Atteints, MAF_Non_atteints)]
  rm(freq)

  # -----
  # Construction du tableau des variants candidats
  # -----
  candidats <- merge(critere_dir[, .(variant_id, probe_id, symbol, directions, n_outliers, n_non_outliers)], 
                outliers_porteurs, by = c("variant_id", "probe_id"), allow.cartesian = TRUE)
  candidats <- merge(candidats, dt_freq, by = "variant_id", all.x = TRUE, )

  # --Filtrage pour les variants rares finaux dans candidats
  fichiers_map <- list.files(filtre_dir, pattern = ".*_chr.*\\.map$", full.names = TRUE)
  variants_rares_finaux <- unique(rbindlist(lapply(fichiers_map, fread, header = FALSE, select = 2))$V2)
  candidats <- candidats[variant_id %in% variants_rares_finaux]

  # --Mise en forme finale
  setnames(candidats, "directions", "direction_expression")
  setorder(candidats, symbol, probe_id, variant_id, IID)

  return(candidats[, .(variant_id, probe_id, symbol, direction_expression, 
                        n_outliers_porteurs = n_outliers, 
                        individu_outlier = IID, 
                        Z_outlier = Z, 
                        MAF_Atteints, MAF_Non_atteints)])

}

# -----------------------------------
# Boucle sur toutes les combinaisons
# -----------------------------------

for (fenetre in fenetres) {
  for (version_z in names(scores_z)) {
    
    data_out <- "/home/chloev/links/projects/def-bureau/chloev/liste_variants/resultats"
    out_file <- file.path(data_out, sprintf("variants_candidats_%s_%s.txt", fenetre, version_z))
    resultat <- filtrer_variants(fenetre, version_z, scores_z[[version_z]])
    
    fwrite(resultat, out_file, sep = "\t", quote = FALSE)

    rm(resultat)
    gc()
  }
}