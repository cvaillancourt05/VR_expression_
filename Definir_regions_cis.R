library(biomaRt)
library(dplyr)
library(readr)

# =======================
# Chargement des données
# =======================

genes <- "/home/chloev/links/projects/def-bureau/expression_genes/positions_SNPs_et_genes/liste_genes/genes_sondes.txt"
genotypes <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/genotyping/omni_grch38.bim"

# ======================================
# Création du fichier BED pour bedtools
# ======================================
all_variants_global <- "all_variants_global.bed"
cmd_bim_to_bed <- sprintf("awk 'print \"chr\"$1\"\t\"($4-1)\"\t\"$4\"\t\"$2}' %s > %s", genotypes, all_variants_global)
system(cmd_bim_to_bed)

# =============================
# Récupération des coordonnées
# =============================

gene_list <- read_lines(genes) %>% trimws() %>% unique()
gene_list <- gene_list[gene_list != ""]

ensembl <- useMart("ensembl", dataset = "hsapiens_genes_ensembl", host = "https://ensembl.org")
queries <- getBM(attributes = c('external_gene_name', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'link_external_gene_name',
  values = gene_list,
  mart = ensembl
)

valid_chromosomes <- c(1:22, "X", "Y")
queries_clean <- queries %>%
  filter(chromosome_name %in% valid_chromosomes) %>%
  mutate(chr = paste0("chr", chromosome_name)) %>%
  rename(gene_name = external_gene_name, start = start_position, end = end_position) %>%
  select(gene_name, chr, start, end)

# =============
# Intersection 
# =============
  
fenetres <- c(10000, 50000) # Fenêtre Mazzarotto et Chagnon
