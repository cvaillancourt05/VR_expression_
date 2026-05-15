# Préparation des fichiers nécessaires pour TensorQTL

library(data.table)
library(readxl)
source("liftOver.R")

# =============================
# Matrice d'expression en hg38 
# =============================

# Fusion avec l'expression
expr_data <- read.delim("data/Adjusted_expression_values.txt", sep = "\t", check.names = FALSE)
# Harmonisation du nom de la première colonnne
colnames(expr_data)[1] <- "probe_id"
# Fusion
expr_hg38 <- merge(lifted, expr_data, by = "probe_id")

# Formattage pour tensorQTL
cols_individu <- setdiff(colnames(expr_hg38), c("probe_id", "illumina_id", "symbol", "chr_hg38", "start_hg38", "end_hg38"))
tensorqtl_bed <- expr_hg38[, c("chr_hg38", "start_hg38", "end_hg38", "probe_id", cols_individu)]
colnames(tensorqtl_bed)[1:4] <- c("#Chr", "start", "end", "ID")
tensorqtl_bed$`#Chr` <- gsub("chr", "", tensorqtl_bed$`#Chr`)
tensorqtl_bed <- tensorqtl_bed[order(tensorqtl_bed$`#Chr`, tensorqtl_bed$start), ]

# Sauvegarde
write.table(tensorqtl_bed, "results/eQTL/expression_hg38.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# ============
# Covariables 
# ============

# Chargement des données du phénotype mis à jour
pheno <- read.table("data/GCbroad_plink.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno) <- c("FID", "IID", "Diagnostic")
pheno$IID <- as.character(pheno$IID)
rownames(pheno) <- pheno$IID

# Éliminer les zéros 
pheno <- pheno[pheno$Diagnostic != 0, ]
# Recodage binaire : Atteint (1), Non-atteint (0)
pheno$statut_Maladie <- ifelse(pheno$Diagnostic == 2, 1, 0)

# Chargement des covariables
df_covars <- read_excel("data/attributs_sujets.xlsx", sheet = 1)
df_covars <- as.data.frame(df_covars)
df_covars$subid <- as.character(df_covars$subid)
rownames(df_covars) <- df_covars$subid

df_covars$Sexe_num <- ifelse(df_covars$sexe == "M", 1, 0)
df_covars$Batch_num <- as.numeric(as.factor(df_covars$batch))

# Intersection
cols_individu_expr <- setdiff(colnames(expr_hg38), c("probe_id", "illumina_id", "symbol", "chr_hg38", "start_hg38", "end_hg38"))

# L'individu doit exister dans tous les fichiers
ind_existants <- intersect(cols_individu_expr, intersect(pheno$IID, df_covars$subid))

# Formatage de la matrice d'expression BED
tensorqtl_bed <- expr_hg38[, c("chr_hg38", "start_hg38", "end_hg38", "probe_id", ind_existants)]
colnames(tensorqtl_bed)[1:4] <- c("#Chr", "start", "end", "ID")
tensorqtl_bed$`#Chr` <- gsub("chr", "", tensorqtl_bed$`#Chr`)
tensorqtl_bed <- tensorqtl_bed[order(tensorqtl_bed$`#Chr`, tensorqtl_bed$start), ]

# Construction de la matrice des covariables
pheno_aligne <- pheno[ind_existants, ]
covars_aligne <- df_covars[ind_existants, ]

df_covariables <- data.frame(
  Statut_Maladie = pheno_aligne$statut_Maladie,
  row.names = ind_existants
)

# Transposition pour TensorQTL
covariates_tensorqtl <- t(df_covariables)

# Sauvegarde
rownames(covariates_tensorqtl) <- c("Phenotype")

# Sauvegarde forcée avec l'identifiant de colonne initial 'id' (Requis par TensorQTL)
write.table(data.frame(id = rownames(covariates_tensorqtl), covariates_tensorqtl, check.names = FALSE), "results/eQTL/covariates.txt", sep = "\t", 
                        row.names = FALSE, col.names = TRUE, quote = FALSE)