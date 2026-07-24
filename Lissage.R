# --------------------------------------------------------------------------------------------
# Lissage.R
# --------------------------------------------------------------------------------------------

library(mgcv)
library(dplyr)

dir_scores <- "results/scores_iogc"
dir_sortie <- "results/asso_pheno"

phenotypes <- c("Pheno_Global", "Pheno_SZ", "Pheno_BP")

# ----------------------------------------------
# Ajustement GAM + graphique de l'effet partiel
# ----------------------------------------------
fonction_gam_lissage <- function(data, pheno, score_col, titre) {

  model_data <- data[!is.na(data[[pheno]]) & !is.na(data[[score_col]]), ]

  form_gam <- as.formula(paste0(pheno, " ~ s(", score_col, ", k = 10)"))
  m <- gam(form_gam, data = model_data, family = binomial(link = "logit"), method = "REML")

  s <- summary(m)
  pv <- s$s.table[1, "p-value"]

  # --seWithMean = TRUE : IC incluant l'incertitude sur la moyenne globale (recommandé, Wood)
  plot(m, shade = TRUE, shade.col = "gray80", rug = TRUE, pch = 19, cex = 0.4,
       seWithMean = TRUE, xlab = "Score IOGC", ylab = "Effet partiel", main = titre)

  data.frame(
    Method = "GAM", Phenotype = pheno, Variable = score_col,
    OR = NA, P_Value = round(pv, 5), CI_Lower = NA, CI_Upper = NA, AUC = NA,
    N_cas = sum(model_data[[pheno]] == 1, na.rm = TRUE), N_total = nrow(model_data), stringsAsFactors = FALSE
  )
}

# ------------------------------------
# Chargement des fichiers phénotypes
# ------------------------------------
lire_pheno <- function(fichier, nom_col) {
  d <- read.table(fichier, header = FALSE, stringsAsFactors = FALSE)
  out <- data.frame(FID = as.character(d$V1), IID = as.character(d$V2),
                     val = ifelse(d$V3 %in% c(1, 2), d$V3 - 1, NA), stringsAsFactors = FALSE)
  names(out)[3] <- nom_col
  out
}

pheno_global <- lire_pheno("data/GCbroad_plink.txt", "Pheno_Global")
pheno_sz <- lire_pheno("data/SZbroad_plink.txt", "Pheno_SZ")
pheno_bp <- lire_pheno("data/BPbroad_plink.txt", "Pheno_BP")

pheno_dt <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID"), all = TRUE),
                    list(pheno_global, pheno_sz, pheno_bp))

# --------------------------
# Boucle sur les conditions
# --------------------------
fichiers_scores <- list.files(dir_scores, pattern = "^score_iogc.*_atteints.*\\.txt$", full.names = TRUE)
resultats <- list()

for (f_path in fichiers_scores) {

  condition <- gsub("\\.txt$", "", sub("score_iogc_", "", basename(f_path)))
  dt_score <- read.delim(f_path, stringsAsFactors = FALSE)
  dt_score$FID <- as.character(dt_score$FID); dt_score$IID <- as.character(dt_score$IID)
  data_clean <- merge(dt_score, pheno_dt, by = c("FID", "IID"))

  # --Un seul fichier PNG par condition, 3 panneaux (1 par phénotype)
  png(file.path(dir_sortie, sprintf("gam_%s.png", condition)), width = 1350, height = 450)
  par(mfrow = c(1, 3))

  for (nom_pheno in phenotypes) {

    n_cas <- sum(data_clean[[nom_pheno]] == 1, na.rm = TRUE)
    if (is.na(n_cas) || n_cas < 10) {
      plot.new(); title(main = paste0(nom_pheno, " : trop peu de cas"))
      next
    }

    res <- tryCatch(fonction_gam_lissage(data_clean, nom_pheno, "score_iogc",
                                          titre = paste0(condition, " - ", nom_pheno)),
                     error = function(e) { message("Échec GAM ", condition, "/", nom_pheno, " : ", e$message); NULL })

    if (!is.null(res)) {
      res$Condition <- condition
      resultats[[paste(condition, nom_pheno, sep = "__")]] <- res
    }
  }

  dev.off()
}

# -----------
# Sauvegarde
# -----------
resultats_finaux <- bind_rows(resultats) %>% select(Condition, Method, Phenotype, P_Value, N_cas, N_total)
write.csv(resultats_finaux, file.path(dir_sortie, "gam_summary.csv"), row.names = FALSE)