# Filtrage des variants candidats

library(data.table)

# ------------------------------------------
# Paramètres - à modifier selon le contexte
# ------------------------------------------

fenetres <- c("10kb", "50kb")
seuil_z <- 2

data_dir <- "/home/chloev/links/projects/def-bureau/chloev/liste_variants/sorties/merge"
scores_z <- c("avec_eQTL_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_avec_eQTL_ref_atteints.txt", 
            "avec_eQTL_non_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_avec_eQTL_ref_non_atteints.txt",
            "sans_eQTL_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_sans_eQTL_ref_atteints.txt",
            "sans_eQTL_non_atteints" = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/Zscores/Z_sans_eQTL_ref_non_atteints.txt")

# -----------------------------------
# Fonction pour filtrer les variants
# -----------------------------------

filtrer_variants <- function(fenetre, version_z, zscore_file) {

  # Intersection variant-sonde
  regions_dt <- fread(file.path(data_dir, sprintf("variants_dans_regions_%s.bed", fenetre)), header = FALSE, col.names = c("variant_id", "probe_gene"))
  regions_dt[, c("probe_id", "symbol") := tstrsplit(probe_gene, "__", fixed = TRUE)]
  regions_dt[, probe_gene := NULL]

  # Matrice de génotypes
  cols_meta <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  raw_files <- list.files(data_dir, pattern = sprintf("genotypes_%s_chr.*\\.raw$", fenetre), full.names = TRUE)
  geno_list <- lapply(raw_files, function(f) {
    dt <- fread(f)
    if (ncol(dt) <= 6) return(NULL)
    dt[, IID := as.character(IID)]
    cols_geno <- setdiff(colnames(dt), cols_meta)
    setnames(dt, cols_geno, sub("_[A-Za-z0-9]+$", "", cols_geno))
    dt[, c("IID", setdiff(colnames(dt), cols_meta)), with = FALSE]
  })
  geno_list <- geno_list[!sapply(geno_list, is.null)]

  geno <- Reduce(function(a, b) merge(a, b, by = "IID", all = TRUE), geno_list)
  geno_long <- melt(geno, id.vars = "IID", variable.name = "variant_id", value.name = "genotype")
  geno_long[, variant_id := as.character(variant_id)]
  geno_porteurs <- geno_long[!is.na(genotype) & genotype > 0]

  # Fréquences alléliques
  freq <- fread(file.path(data_dir, sprintf("freq_variants_%s.afreq", fenetre)))
  setnames(freq, old = c("CHROM", "ID"), new = c("CHROM", "variant_id"), skip_absent= TRUE)
  dt_freq <- freq[, .(variant_id, REF, ALT, ALT_FREQ)]

  # Scores Z
  zscores <- fread(zscore_file)
  setnames(zscores, 1, "probe_id")
  z_long <- melt(zscores, id.vars = "probe_id", variable.name = "IID", value.name = "Z")
  z_long[, IID := as.character(IID)]
  z_long[, is_outlier := abs(Z) >= seuil_z]
  z_long[, direction := fifelse(Z >= seuil_z, "UP", fifelse(Z <= -seuil_z, "DOWN", "none"))]

  # Table de travail
  work <- merge(regions_dt, geno_porteurs, by = "variant_id")
  work <- merge(work, z_long[, .(probe_id, IID, Z, is_outlier, direction)], by = c("probe_id", "IID"), all.x = TRUE)

  # Critère de présence
  resume <- work[!is.na(is_outlier), .(n_outliers = sum(is_outlier), n_non_outliers = sum(!is_outlier), directions = paste(unique(direction[is_outlier]), collapse = ",")), by = .(variant_id, probe_id, symbol)]
  critere_pres <- resume[n_outliers >= 1 & n_non_outliers == 0]

  # Critère de direction cohérente
  critere_dir <- critere_pres[directions %in% c("UP", "DOWN")]

  work[, combo_key := paste(variant_id, probe_id)]
  critere_dir[, combo_key := paste(variant_id, probe_id)]

  outliers_porteurs <- work[is_outlier == TRUE & combo_key %in% critere_dir$combo_key, .(individus_outliers = paste(IID, collapse = ";")), by = .(variant_id, probe_id)]

  z_outliers <- work[is_outlier == TRUE & combo_key %in% critere_dir$combo_key, .(Z_outliers = paste(round(Z, 3), collapse = ";")), by = .(variant_id, probe_id)]

  # Construction du tableau final des variants candidats
  candidats <- merge(critere_dir, outliers_porteurs, by = c("variant_id", "probe_id"))
  candidats <- merge(candidats, z_outliers, by = c("variant_id", "probe_id"))
  candidats <- merge(candidats, dt_freq, by = "variant_id", all.x = TRUE)

  setnames(candidats, "directions", "direction_expression")
  setorder(candidats, symbol, probe_id, variant_id)

  return(candidats[, .(variant_id, probe_id, symbol, direction_expression, n_outliers_porteurs = n_outliers, individus_outliers, Z_outliers, REF, ALT, ALT_FREQ)])

}

# -----------------------------------
# Boucle sur toutes les combinaisons
# -----------------------------------

for (fenetre in fenetres) {
  for (version_z in names(scores_z)) {
    
    out_file <- file.path(data_dir, sprintf("variants_candidats_%s_%s.txt", fenetre, version_z))
    resultat <- filtrer_variants(fenetre, version_z, scores_z[[version_z]])
    
    fwrite(resultat, out_file, sep = "\t", quote = FALSE

    )
  }
}