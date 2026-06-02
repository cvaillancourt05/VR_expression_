# ---------------------------------------------------------------------
# prep_tensorQTL.R
# Préparation des fichiers d'entrée pour TensorQTL (analyse cis-eQTL)
#
# Entrées : 
#  - Fichier de la matrice d'expression ajustée
#  - Fichier de génotypes
# Sorties : 
#  - Vecteur de statut diagnostic -- phenotypes.txt
#  - Fichier de la matrice d'expression au format BED -- expression_hg38.bed.gz
#  - Index Tabix du fichier BED -- expression_hg38.bed.gz.tbi
# ---------------------------------------------------------------------

library(data.table)
library(readxl)
library(Rsamtools)
# --Pour charger la table 'lifted'
source("liftOver.R")

# --------------------------------------------------
# Fusion de l'expression et les coordonnées en hg38 
# --------------------------------------------------

expr_data <- read.delim("data/Adjusted_expression_values.txt", sep = "\t", check.names = FALSE)
colnames(expr_data)[1] <- "probe_id" # --Harmonisation du nom de la première colonnne

expr_hg38 <- merge(lifted, expr_data, by = "probe_id")

# ----------------------------------------------------------
# Chargement du phénotype et harmonisation des identifiants  
# ----------------------------------------------------------

pheno <- read.table("data/GCbroad_plink.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno) <- c("FID", "IID", "Diagnostic")
pheno$IID <- as.character(pheno$IID)
pheno$FID <- as.character(pheno$FID)

# --Création de l'ID unique au format FID_IID 
pheno$ID <- paste(pheno$FID, pheno$IID, sep = "_")
rownames(pheno) <- pheno$ID

# --Éliminer les individus au status inconnu (Diagnostic == 0)
pheno <- pheno[pheno$Diagnostic != 0, ]
# --Recodage binaire : Atteint (1), Non-atteint (0)
pheno$statut_Maladie <- ifelse(pheno$Diagnostic == 2, 1, 0)

# --Harmonisation des ID
cols_base <- c("probe_id", "illumina_id", "symbol", "chr_hg38", "start_hg38", "end_hg38")
ids_expression <- setdiff(colnames(expr_hg38), cols_base)
map_ids <- setNames(pheno$ID, pheno$IID)
new_ids <- map_ids[ids_expression]

ids_valide <- ids_expression[!is.na(new_ids)]
new_ids_valides <- new_ids[!is.na(new_ids)]

colnames(expr_hg38)[match(ids_valide, colnames(expr_hg38))] <- new_ids_valides

# --Intersection
# Garder uniquement les individus présents dans l'expression et le phénotype
ids_expression <- setdiff(colnames(expr_hg38), cols_base)
ind_existants <- intersect(ids_expression, pheno$ID)

# ---------------------------------------------------
# Formatage de la matrice d'expression au format BED
# ---------------------------------------------------

tensorqtl_bed <- expr_hg38[, c("chr_hg38", "start_hg38", "end_hg38", "probe_id", ind_existants)]
colnames(tensorqtl_bed)[1:4] <- c("#Chr", "start", "end", "ID")

# --TensorQTL attend des numéros de chromosome sans préfie 'chr'
tensorqtl_bed$`#Chr` <- gsub("chr", "", tensorqtl_bed$`#Chr`)
tensorqtl_bed$start <- as.numeric(tensorqtl_bed$start)
tensorqtl_bed$end <- as.numeric(tensorqtl_bed$end)

# --Tri par chromosome puis position
tensorqtl_bed$chr_num <- suppressWarnings(as.numeric(tensorqtl_bed$`#Chr`))
tensorqtl_bed$chr_num[tensorqtl_bed$`#Chr` == "X"] <- 23
tensorqtl_bed$chr_num[tensorqtl_bed$`#Chr` == "Y"] <- 24
tensorqtl_bed$chr_num[tensorqtl_bed$`#Chr` == "M" | tensorqtl_bed$`#Chr` == "MT"] <- 25
tensorqtl_bed$chr_num[is.na(tensorqtl_bed$chr_num)] <- 99 # --chromosomes non standard en dernier

tensorqtl_bed <- tensorqtl_bed[order(tensorqtl_bed$chr_num, tensorqtl_bed$start, tensorqtl_bed$end), ]
tensorqtl_bed$chr_num <- NULL

# -------------------------------------------------------
# Construction de la matrice des covariables (phénotype)
# -------------------------------------------------------

# --Aligner le phénotype sur les individus retenus dans la matrice d'expression
pheno_aligne <- pheno[ind_existants, ]

# --Transposer 
pheno_t <- data.frame(matrix(pheno_aligne$statut_Maladie, nrow = 1))
colnames(pheno_t) <- ind_existants
rownames(pheno_t) <- "Phenotype"

pheno_tensorqtl <- cbind(`#id` = rownames(pheno_t), pheno_t)

# ----------------------------------
# Sauvegarde des fichiers de sortie
# ----------------------------------

# --Fichier de phénotypes
write.table(pheno_tensorqtl, "results/phenotypes.txt", sep = "\t", 
                        row.names = FALSE, col.names = TRUE, quote = FALSE)

# --Fichier BED brut intermédiaire
write.table(tensorqtl_bed, "results/expression_hg38.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# --Compression en BGZF + indexation Tabix
out <- bgzip("results/expression_hg38.bed", dest = "results/expression_hg38.bed.gz", overwrite = TRUE)
indexTabix(out, format = "bed")