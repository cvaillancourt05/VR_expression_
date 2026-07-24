# --------------------------------------------------------------------------------------------
# Association_phenotypique_p90.R
# Adaptation de la méthodologie de Mazzarotto et al. 2025 à la cohorte Chagnon
# Analyse du fardeau de variants rares IOGC dichotomisé au 90e percentile (P90 vs reste)
# Pas de PRS ni de prévalence populationnelle (K) : métriques basées sur l'échelle observée.
#
# Entrées :
#  - results/scores_iogc/score_iogc_*_atteints*.txt
#  - data/GCbroad_plink.txt, SZbroad_plink.txt, BPbroad_plink.txt
#
# Sortie :
#  - results/asso_pheno/associations_p90.csv
# --------------------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(geepack)
library(readxl)
library(data.table)


dir_scores <- "results/scores_iogc"
dir_sortie <- "results/asso_pheno"

# --------------------------------------------------------------------
# Dichotomisation du score au 90e percentile (P90)
# --------------------------------------------------------------------
dichotomiser_p90 <- function(data, score_col) {
  seuil <- quantile(data[[score_col]], 0.90, na.rm = TRUE)
  data$groupe_iogc <- factor(ifelse(data[[score_col]] > seuil, "haut", "bas"), levels = c("bas", "haut"))
  attr(data, "seuil_p90") <- seuil
  data
}

# --------------------------------------------------------------------
# Fonction de régression logistique standard (GLM)
# --------------------------------------------------------------------
fonction_glm <- function(data, pheno, score_col) {
  model_data <- data %>% 
    drop_na(all_of(c(pheno, score_col, "Sexe"))) %>% 
    dichotomiser_p90(score_col)
  
  model <- glm(as.formula(paste0(pheno, " ~ groupe_iogc + Sexe")), data = model_data, family = binomial(link = "logit"))
  summaryModel <- summary(model)
  coef_df <- as.data.frame(summaryModel$coefficients)
  
  if(!"groupe_iogchaut" %in% rownames(coef_df)) return(NULL)
  coef_row <- coef_df["groupe_iogchaut", ]

  data.frame(
    Method = "GLM", Phenotype = pheno, Variable = "P90_vs_reste",
    OR = round(exp(coef_row[["Estimate"]]), 3),
    P_Value = round(coef_row[["Pr(>|z|)"]], 5),
    CI_Lower = round(exp(coef_row[["Estimate"]] - 1.96 * coef_row[["Std. Error"]]), 3),
    CI_Upper = round(exp(coef_row[["Estimate"]] + 1.96 * coef_row[["Std. Error"]]), 3),
    N_cas = sum(model_data[[pheno]] == 1, na.rm = TRUE), N_total = nrow(model_data),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------
# Fonction des équations d'estimation généralisées (GEE) - Structure familiale
# --------------------------------------------------------------------
fonction_gee <- function(data, pheno, score_col) {
  model_data <- data %>% 
    drop_na(all_of(c(pheno, "FID", score_col, "Sexe"))) %>% 
    dichotomiser_p90(score_col) %>% 
    arrange(FID)
  
  model <- geeglm(as.formula(paste0(pheno, " ~ groupe_iogc + Sexe")), data = model_data, id = factor(FID), 
                  family = binomial(link = "logit"), corstr = "independence")
  summaryModel <- summary(model)
  coef_df <- as.data.frame(summaryModel$coefficients)
  
  if(!"groupe_iogchaut" %in% rownames(coef_df)) return(NULL)
  coef_row <- coef_df["groupe_iogchaut", ]
  
  a <- coef_row[["Estimate"]]
  b <- coef_row[["Std.err"]]

  data.frame(
    Method = "GEE", Phenotype = pheno, Variable = "P90_vs_reste", 
    OR = round(exp(a), 3),
    P_Value = round(coef_row[["Pr(>|W|)"]], 5),
    CI_Lower = round(exp(a - 1.96 * b), 3),
    CI_Upper = round(exp(a + 1.96 * b), 3),
    N_cas = sum(model_data[[pheno]] == 1, na.rm = TRUE), N_total = nrow(model_data),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------
# Chargement et Standardisation des fichiers phénotypes
# --------------------------------------------------------------------
pheno_global <- read_table("data/GCbroad_plink.txt", show_col_types = FALSE, col_names = FALSE) %>%
  mutate(FID = as.character(X1), IID = as.character(X2), Pheno_Global = if_else(X3 %in% c(1, 2), X3 - 1, NA_real_)) %>%
  select(FID, IID, Pheno_Global)

pheno_sz <- read_table("data/SZbroad_plink.txt", show_col_types = FALSE, col_names = FALSE) %>%
  mutate(FID = as.character(X1), IID = as.character(X2), Pheno_SZ = if_else(X3 %in% c(1, 2), X3 - 1, NA_real_)) %>%
  select(FID, IID, Pheno_SZ)

pheno_bp <- read_table("data/BPbroad_plink.txt", show_col_types = FALSE, col_names = FALSE) %>%
  mutate(FID = as.character(X1), IID = as.character(X2), Pheno_BP = if_else(X3 %in% c(1, 2), X3 - 1, NA_real_)) %>%
  select(FID, IID, Pheno_BP)

pheno_all <- pheno_global %>%
  full_join(pheno_sz, by = c("FID", "IID")) %>%
  full_join(pheno_bp, by = c("FID", "IID"))

covar_sexe <- as.data.table(read_excel("data/attributs_sujets.xlsx")) %>%
  mutate(IID = as.character(subid), Sexe = as.factor(sexe.x)) %>%
  select(IID, Sexe)
pheno_all <- pheno_all %>% left_join(covar_sexe, by = "IID")

phenotypes <- c("Pheno_Global", "Pheno_SZ", "Pheno_BP")

# --------------------------------------------------------------------
# Boucle principale d'analyse pour chaque condition
# --------------------------------------------------------------------

fichiers_scores <- list.files(dir_scores, pattern = "^score_iogc.*_atteints.*\\.txt$", full.names = TRUE)
resultats <- list()

for (f_path in fichiers_scores) {

  condition <- gsub("\\.txt$", "", sub("score_iogc_", "", basename(f_path)))
  dt_score <- read.delim(f_path, stringsAsFactors = FALSE)
  dt_score$FID <- as.character(dt_score$FID)
  dt_score$IID <- as.character(dt_score$IID)
  
  data_clean <- merge(dt_score, pheno_all, by = c("FID", "IID"))

  for (nom_pheno in phenotypes) {
    res_glm <- tryCatch(fonction_glm(data_clean, nom_pheno, "score_iogc"), 
                        error = function(e) { message("Échec GLM ", condition, "/", nom_pheno, " : ", e$message); NULL})
    res_gee <- tryCatch(fonction_gee(data_clean, nom_pheno, "score_iogc"),
                        error = function(e) { message("Échec GEE ", condition, "/", nom_pheno, " : ", e$message); NULL })
  
    res_stock <- do.call(rbind, list(res_glm, res_gee))
    if (!is.null(res_stock) && nrow(res_stock) > 0) {
      res_stock$Condition <- condition
      resultats[[paste(condition, nom_pheno, sep = "__")]] <- res_stock
    }
  }
}

# --------------------------------------------------------------------
# Sauvegarde des résultats
# --------------------------------------------------------------------

if (length(resultats) > 0) {
  resultats_finaux <- bind_rows(resultats) %>%
  select(Condition, Method, Phenotype, OR, P_Value, CI_Lower, CI_Upper, N_cas, N_total)
  write_csv(resultats_finaux, file.path(dir_sortie, "associations_p90_covar_sexe.csv"))
}