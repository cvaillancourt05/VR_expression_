# Préparation des fichiers nécessaires pour TensorQTL

library(data.table)
library(readxl)
library(Rsamtools)
source("liftOver.R")

# =============================
# Matrice d'expression en hg38 
# =============================

# Fusion avec l'expression
expr_data <- read.delim("/home/chloev/links/projects/def-bureau/expression_genes/Adjusted_expression_values.txt", sep = "\t", check.names = FALSE)
# Harmonisation du nom de la première colonnne
colnames(expr_data)[1] <- "probe_id"

sondes_expr <- read.table("/home/chloev/links/projects/def-bureau/expression_genes/liste_sondes_expr_GC.txt", header = FALSE, stringsAsFactors = FALSE)$V1
expr_data <- expr_data[expr_data$probe_id %in% sondes_expr, ]

# Fusion
expr_hg38 <- merge(lifted, expr_data, by = "probe_id")

# ===========================
# Intersection et alignement  
# ===========================

# Chargement des données du phénotype mis à jour
pheno <- read.table("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/GCbroad_plink.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pheno) <- c("FID", "IID", "Diagnostic")
pheno$IID <- as.character(pheno$IID)
rownames(pheno) <- pheno$IID

# Éliminer les zéros 
pheno <- pheno[pheno$Diagnostic != 0, ]
# Recodage binaire : Atteint (1), Non-atteint (0)
pheno$statut_Maladie <- ifelse(pheno$Diagnostic == 2, 1, 0)

# Intersection
cols_individu_expr <- setdiff(colnames(expr_hg38), c("probe_id", "illumina_id", "symbol", "chr_hg38", "start_hg38", "end_hg38"))

# L'individu doit exister dans tous les fichiers
ind_existants <- intersect(cols_individu_expr, pheno$IID)

plink_fam <- read.table("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/genotyping/omni_grch38.fam", header = FALSE, stringsAsFactors = FALSE)
plink_ids <- as.character(plink_fam$V2)

ind_existants <- intersect(intersect(cols_individu_expr, pheno$IID), "\n")

# =========================================
# Formatage de la matrice d'expression BED
# =========================================

tensorqtl_bed <- expr_hg38[, c("chr_hg38", "start_hg38", "end_hg38", "probe_id", ind_existants)]
colnames(tensorqtl_bed)[1:4] <- c("#Chr", "start", "end", "ID")
tensorqtl_bed$`#Chr` <- gsub("chr", "", tensorqtl_bed$`#Chr`)

tensorqtl_bed$start <- as.numeric(tensorqtl_bed$start)
tensorqtl_bed$end <- as.numeric(tensorqtl_bed$end)

tensorqtl_bed$chr_num <- suppressWarnings(as.numeric(tensorqtl_bed$`#Chr`))
tensorqtl_bed$chr_num[tensorqtl_bed$`#Chr` == "X"] <- 23
tensorqtl_bed$chr_num[tensorqtl_bed$`#Chr` == "Y"] <- 24
tensorqtl_bed$chr_num[tensorqtl_bed$`#Chr` == "M" | tensorqtl_bed$`#Chr` == "MT"] <- 25
tensorqtl_bed$chr_num[is.na(tensorqtl_bed$chr_num)] <- 99

tensorqtl_bed <- tensorqtl_bed[order(tensorqtl_bed$chr_num, tensorqtl_bed$start, tensorqtl_bed$end), ]
tensorqtl_bed$chr_num <- NULL


# ==========================================
# Construction de la matrice des covariables
# ==========================================

pheno_aligne <- pheno[ind_existants, ]

pheno_t <- data.frame(matrix(pheno_aligne$statut_Maladie, nrow = 1))
colnames(pheno_t) <- ind_existants
rownames(pheno_t) <- "Phenotype"

pheno_tensorqtl <- cbind(`#id` = rownames(pheno_t), pheno_t)

# ===========
# Sauvegarde
# ===========

# Sauvegarde du fichier de phénotypes
write.table(pheno_tensorqtl, "/home/chloev/links/projects/def-bureau/expression_genes/resultats/phenotypes.txt", sep = "\t", 
                        row.names = FALSE, col.names = TRUE, quote = FALSE)

# Sauvegarde du fichier BED brut
write.table(tensorqtl_bed, "/home/chloev/links/projects/def-bureau/expression_genes/resultats/expression_hg38.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Sauvegarde en .bed.gz 
#tensorqtl_bed_zip <- bgzip("results/eQTL/expression_hg38.bed", dest= "results/eQTL/expression_hg38.bed.gz", overwrite = TRUE)
#indexTabix(tensorqtl_bed_zip, format = "bed")

out <- bgzip("/home/chloev/links/projects/def-bureau/expression_genes/resultats/expression_hg38.bed", dest = "/home/chloev/links/projects/def-bureau/expression_genes/resultats/expression_hg38.bed.gz", overwrite = TRUE)
indexTabix(out, format = "bed")