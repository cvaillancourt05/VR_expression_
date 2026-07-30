# ---------------------------------------------------------------------------------
# filter_samples.R
# Filtrage de expression_hg38.bed.gz et de phenotypes.txt 
# ---------------------------------------------------------------------------------


base_dir <- "/lustre09/project/6000443/expression_genes/resultats/retrait_effet_eqtl"

# --Liste des échantillons du génotype final
fam <- read.table(file.path(base_dir, "genotypes.fam"),
                   stringsAsFactors = FALSE)
samples_geno <- trimws(as.character(fam$V2))

# --Fonction qui lit le header d'un fichier
read_header_clean <- function(path, is_gz = FALSE) {
    con <- if (is_gz) gzfile(path) else file(path)
    header_line <- readLines(con, n = 1)
    close(con)
    header_line <- gsub("\r$", "", header_line)    
    header_line <- gsub("^\ufeff", "", header_line) 
    trimws(strsplit(header_line, "\t")[[1]])
}

# --Filtrage de phenotypes.txt (covariates x samples)
header_pheno <- read_header_clean(file.path(base_dir, "phenotypes.txt"))
sample_cols_pheno <- header_pheno[-1]  # toutes les colonnes sauf #id

covariates_df <- read.table(file.path(base_dir, "phenotypes.txt"),
                             sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
                             comment.char = "")
colnames(covariates_df) <- sample_cols_pheno 

samples_communs <- intersect(trimws(colnames(covariates_df)), samples_geno)

covariates_df <- covariates_df[, samples_communs]
write.table(covariates_df, file.path(base_dir, "phenotypes_filtered.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# --Filtrage de expression_hg38.bed.gz (chr, start, end, id + samples)
header_bed <- read_header_clean(file.path(base_dir, "expression_hg38.bed.gz"), is_gz = TRUE)

pheno_bed <- read.table(gzfile(file.path(base_dir, "expression_hg38.bed.gz")),
                         sep = "\t", header = TRUE, check.names = FALSE, comment.char = "")
colnames(pheno_bed) <- header_bed  

cols_a_garder <- c(colnames(pheno_bed)[1:4], samples_communs)
pheno_bed <- pheno_bed[, cols_a_garder]

# -- Écriture en texte brut (non compressé) ; bgzip + tabix se fait après, en bash
write.table(pheno_bed, file.path(base_dir, "expression_hg38_filtered.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE)