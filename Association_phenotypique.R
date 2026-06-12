# --------------------------------------------------------------------------------------------
# Association_phenotypique.R
# Réalise les associations statistiques entre les scores IOGC et les phénotypes (SZ, BD, CL)
#
# Entrées :
#  - score_iogc.csv
#
# Sortie :
#  - asssociations_finales_iogc.csv
# -------------------------------------------------------------------------------------------

library(geepack)
library(pROC)
library(dplyr)
library(tidyr)
library(readr)

# ------------------------------------------------
# Fonction de régression logistique standard (GLM)
# ------------------------------------------------

fonction_glm <- function(data, pheno, score_col) {

  variables <- score_col
  form_str <- paste0(pheno, " ~ ", paste(variables, collapse = " + "))
  model_data <- data %>% drop_na(all_of(c(pheno, variables)))

  model <- glm(as.formula(form_str), data = model_data, family = binomial(link = "logit"))
  summaryModel <- summary(model)

  if(!score_col %in% rownames(summaryModel$coefficients)) return(NULL)
  coef_row <- summaryModel$coefficients[score_col, ]

  roc_obj <- roc(response = model_data[[pheno]], 
                  predictor = model_data[[score_col]], 
                  direction = "<", levels = c(0, 1), quiet = TRUE)
  
  data.frame(
    Method = "GLM", Phenotype = pheno, Variable = score_col,
    OR = round(exp(coef_row["Estimate"]), 3),
    P_Value = round(coef_row["Pr(>|z|)"], 5),
    CI_Lower = round(exp(coef_row["Estimate"] - 1.96 * coef_row["Std. Error"]), 3), 
    CI_Upper = round(exp(coef_row["Estimate"] + 1.96 * coef_row["Std. Error"]), 3),
    AUC = round(as.numeric(pROC::auc(roc_obj)), 3),
    N_cas = sum(model_data[[pheno]] == 1, na.rm = TRUE), N_total = nrow(model_data),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------
# Fonction des équations d'estimation généralisées (GEE)
# ------------------------------------------------------

fonction_gee <- function(data, pheno, score_col) {

  variables <- score_col
  form_str <- paste0(pheno, " ~ ", paste(variables, collapse = " + "))
  model_data <- data %>% drop_na(all_of(c(pheno, "FID", variables))) %>% arrange(FID)

  model <- geeglm(as.formula(form_str), data = model_data, id = factor(FID), 
            family = binomial(link = "logit"), corstr = "independence", scale.fix = TRUE, scale.value = rep(1, length(model_data$FID)))
  summaryModel <- summary(model)

  if (!score_col %in% rownames(summaryModel$coefficients)) return (NULL)
  coef_row <- summaryModel$coefficients[score_col, ]

  roc_obj <- roc(response = model_data[[pheno]], 
                  predictor = model_data[[score_col]], direction = "<", 
                  levels = c(0, 1), quiet = TRUE)
  
  data.frame(
    Method = "GEE", Phenotype = pheno, Variable = score_col,
    OR = round(exp(coef_row["Estimate"]), 3),
    P_Value = round(coef_row["Pr(>|W|)"], 5),
    CI_Lower = round(exp(coef_row["Estimate"] - 1.96 * coef_row["Std.err"]), 3),
    CI_Upper = round(exp(coef_row["Estimate"] + 1.96 * coef_row["Std.err"]), 3),
    AUC = round(as.numeric(pROC::auc(roc_obj)), 3),
    N_cas = sum(model_data[[pheno]] == 1, na.rm = TRUE), N_total = nrow(model_data),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------
# Chargement des fichiers phénotypes
# -----------------------------------

pheno_global <- read_table("results/phenotypes_formate.txt", show_col_types = FALSE) %>%
  mutate(FID = as.character(FID), IID = as.character(IID), Pheno_Global = as.numeric(Phenotype)) %>%
  select(FID, IID, Pheno_Global)

pheno_sz <- read_table("data/SZbroad_plink.txt", show_col_types = FALSE) %>%
  mutate(FID = as.character(FID), IID = as.character(IID), Pheno_SZ = as.numeric(Phenotype)) %>%
  select(FID, IID, Pheno_SZ)

pheno_bp <- read_table("data/BPbroad_plink.txt", show_col_types = FALSE) %>%
  mutate(FID = as.character(FID), IID = as.character(IID), Pheno_BP = as.numeric(Phenotype)) %>%
  select(FID, IID, Pheno_BP)

# ------------------------------------
# Boucle pour toutes les combinaisons
# ------------------------------------

scores_dir <- "/home/chloev/links/projects/def-bureau/chloev/"
fichiers_scores <- list.files(scores_dir, pattern = "^score_iogc_.*\\.txt$", full.names = TRUE)

phenotypes <- list("Pheno_Global", "Pheno_SZ", "Pheno_BP")
pheno_dt <- reduce(list(pheno_global, pheno_sz, pheno_bp), full_join)

resultats <- list()

for (f_path in fichiers_scores) {

  # --Extraire le nom du fichier
  f_name <- basename(f_path)
  condition <- sub("score_iogc_", "", f_name)
  condition <- gsub("\\.txt$", "", condition)

  # --Charger le fichier
  dt_score <- read_tsv(f_path, show_col_types = FALSE) %>% mutate(FID = as.character(FID), IID = as.character(IID))

  for (nom_pheno in phenotypes) {

    data_clean <- dt_score %>% inner_join(pheno_dt, by = c("FID", "IID"))
    cle_combinaison <- paste(condition, nom_pheno, sep = "__")

    # --GLM
    res_glm <- tryCatch(fonction_glm(data_clean, nom_pheno, "score_iogc"),
                        error = function(e) { message(" Échec GLM pour ", cle_combinaison, " : ", e$message); NULL })
    # --GEE
    res_gee <- tryCatch(fonction_gee(data_clean, nom_pheno, "score_iogc"),
                        error = function(e) { message(" Échec GEE pour ", cle_combinaison, " : ", e$message); NULL })

    # --Stockage des résultats
    res_stock <- bind_rows(res_glm, res_gee)
    if (!is.null(res_stock) && nrow(res_stock) > 0) {
      res_stock$Condition <- condition
      resultats[[cle_combinaison]] <- res_stock
    }
  }
}

# -----------
# Sauvegarde
# -----------
resultats_finaux <- bind_rows(resultats) %>%
  select(Condition, Method, Phenotype, OR, P_Value, CI_Lower, CI_Upper, AUC, N_cas, N_total)

write_csv(resultats_finaux, "results/associations_finales_iogc.csv")