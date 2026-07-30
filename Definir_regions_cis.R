# ------------------------------------------------------------------------------------------------------
# Definir_regions_cis.R
# Définit les régions cis autour des gènes des sondes d'expression et génère les fichiers au format BED
# 
# Entrées :
#  - Annotation des sondes après le liftOver hg19 -> hg38
#  - Liste des sondes retenues après filtrage
#
# Sorties :
#  - Table de correspondance sonde-gène avec coordonnées en hg38 -- probe_gene_map.txt
#  - Régions cis +/- 10kb autour de chaque gène -- region_cis_10kb.bed
#  - Région cis +/- 50kb autour de chaque gène -- region_cis_50kb.bed
# ------------------------------------------------------------------------------------------------------

library(biomaRt)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# -----------------------
# Chargement des données
# -----------------------

lifted <- fread("/lustre09/project/6000443/expression_genes/resultats/HumanHT12v4_hg38_annotations.csv")
sondes_exprimees <- fread("/lustre09/project/6000443/expression_genes/resultats/liste_sondes_exprimees.txt", header = FALSE, col.names = "probe_id")$probe_id

lifted <- lifted[probe_id %in% sondes_exprimees]
lifted <- lifted[!is.na(symbol) & symbol != "" & symbol != "---"]


# ----------------------------------------
# Récupération des coordonnées génomiques
# ----------------------------------------

symboles_uniques <- unique(lifted$symbol)

# --Coordonnées en hg38 depuis TxDb
genes <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
# --Correspondance des symboles
mapping <- select(org.Hs.eg.db, keys = symboles_uniques, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")

mapping <- na.omit(mapping)
mapping <- mapping[!duplicated(mapping$SYMBOL), ]

# --Filtrage sur les gènes présents dans TxDb et conversion en data table
ids_valide <- intersect(mapping$ENTREZID, names(genes))
genes_filtrees <- genes[ids_valide]
gene_coords <- as.data.table(genes_filtrees)
gene_coords$hgnc_symbol <- mapping$SYMBOL[match(gene_coords$gene_id, mapping$ENTREZID)]

# --Mise en forme des coordonnées et filtrage sur les chromosomes standards
gene_coords[, chromosome_name := gsub("chr", "", seqnames)]
gene_coords[, start_position := start]
gene_coords[, end_position := end]
gene_coords <- gene_coords[chromosome_name %in% c(as.character(1:22), "X", "Y")]
gene_coords[, chr := paste0("chr", chromosome_name)]

# --En cas de gène dupliqué, conserver la région la plus longue 
gene_coords[, gene_length := end_position - start_position]
setorder(gene_coords, hgnc_symbol, -gene_length)
gene_coords <- gene_coords[!duplicated(hgnc_symbol) & !is.na(hgnc_symbol)]


# --Table de correspondance sonde-gène avec des coordonnées hg38
probe_gene <- merge(
  lifted[, .(probe_id, illumina_id, symbol)], 
  gene_coords[, .(hgnc_symbol, chr, start_position, end_position)], 
  by.x = "symbol", by.y = "hgnc_symbol", all.x = FALSE)

fwrite(probe_gene, "/lustre09/project/6000443/expression_genes/resultats/probe_gene_map.txt", sep = "\t", quote = FALSE)


# --------------------------------------
# Création du fichier BED pour bedtools
# --------------------------------------

# --Format BED : chr | start (0-based) | end | ID (sonde)
creer_bed <- function(dt, fenetre_kb) {
  fenetre <- fenetre_kb * 1000

  bed <- data.table(
    chr = dt$chr,
    start = pmax(0, dt$start_position - fenetre),
    end = dt$end_position + fenetre,
    ID = paste0(dt$probe_id, "__", dt$symbol)
  )

  # --Tri chromosomique numérique
  bed[, chr_num := suppressWarnings(as.numeric(gsub("chr", "", chr)))]
  bed[chr == "chrX", chr_num := 23]
  bed[chr == "chrY", chr_num := 24]
  setorder(bed, chr_num, start)
  bed[, chr_num := NULL]

  return(bed)
}

# ----------------------------
# Génération des fichiers BED
# ----------------------------

bed_10kb <- creer_bed(probe_gene, fenetre_kb = 10)
bed_50kb <- creer_bed(probe_gene, fenetre_kb = 50)

fwrite(bed_10kb, "/lustre09/project/6000443/expression_genes/resultats/regions_cis_10kb.bed", sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(bed_50kb, "/lustre09/project/6000443/expression_genes/resultats/regions_cis_50kb.bed", sep = "\t", col.names = FALSE, quote = FALSE)