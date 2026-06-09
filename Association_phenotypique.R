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

fonction_glm <- function(data, pheno, score_col, covars) {

  variables <- c(score_col, covars)
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

fonction_gee <- function(data, pheno, score_col, covars) {

  variables <- c(score_col, covars)
  form_str <- paste0(pheno, " ~ ", paste(variables, collapse = " + "))
  model_data <- data %>% drop_na(all_of(c(pheno, "FID", variables))) %>% arrange(FID)

  model <- geeglm(as.formula(form_str), data = model_data, id = factor(FID), 
            family = binomial(link = "logit", corstr = "independence", scale.fix = TRUE, scale.value = 1))
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

# ------------------------------------
# Boucle pour toutes les combinaisons
# ------------------------------------

data <- read_csv("scores_iogc.csv", show_col_types = FALSE)
covariables <- NULL
phenotype = "Phenotype"

col_scores <- grep("IOGC_adj_", colnames(data), value = TRUE)

resultats <- list()
for (score in col_scores) {
  # --GLM
  res_glm <- tryCatch(fonction_glm(data, pheno, score, covariables),
                        error = function(e) { message (" Échec GLM : ", e$message); NULL})
  # --GEE
  res_gee <- tryCatch(fonction_gee(data, pheno, score, covariables),
                        error = function(e) { message (" Échec GEE : ", e$message); NULL})

  resultats[[score]] <- bind_rows(res_glm, res_gee)
}

resultats_finaux <- bind_rows(resultats) %>% 
  mutate(Condition = gsub("IOGC_adj_", "", Variable)) %>% 
  select(Condition, Method, Phenotype, OR, P_Value, CI_Lower, CI_Upper, AUC, N_cas, N_total)

# --Sauvegarde
write_csv(resultats_finaux, "associations_finales_iogc.csv")