library(biomaRt)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# =======================
# Chargement des données
# =======================

lifted <- fread("results/HumanHT12v4_hg38_annotations.csv")
sondes_exprimees <- fread("results/liste_sondes_exprimees.txt", header = FALSE, col.names = "probe_id")$probe_id

lifted <- lifted[probe_id %in% sondes_exprimees]
lifted <- lifted[!is.na(symbol) & symbol != "" & symbol != "---"]


# =============================
# Récupération des coordonnées
# =============================

symboles_uniques <- unique(lifted$symbol)

genes <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
# Correspondance entre ID et symbole
mapping <- select(org.Hs.eg.db, keys = symboles_uniques, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")

mapping <- na.omit(mapping)
mapping <- mapping[!duplicated(mapping$SYMBOL), ]

ids_valide <- intersect(mapping$ENTREZID, names(genes))
genes_filtrees <- genes[ids_valide]

gene_coords <- as.data.table(genes_filtrees)
gene_coords$hgnc_symbol <- mapping$SYMBOL[match(gene_coords$gene_id, mapping$ENTREZID)]

gene_coords[, chromosome_name := gsub("chr", "", seqnames)]
gene_coords[, start_position := start]
gene_coords[, end_position := end]

# Filtrage
gene_coords <- gene_coords[chromosome_name %in% c(as.character(1:22), "X", "Y")]
gene_coords[, chr := paste0("chr", chromosome_name)]

# En cas de gène dupliqué, garder la plus longue région
gene_coords[, gene_length := end_position - start_position]
setorder(gene_coords, hgnc_symbol, -gene_length)
gene_coords <- gene_coords[!duplicated(hgnc_symbol) & !is.na(hgnc_symbol)]


# Pour les sondes
probe_gene <- merge(
  lifted[, .(probe_id, illumina_id, symbol)], 
  gene_coords[, .(hgnc_symbol, chr, start_position, end_position)], 
  by.x = "symbol", by.y = "hgnc_symbol", all.x = FALSE)

fwrite(probe_gene, "results/region_cis/probe_gene_map.txt", sep = "\t", quote = FALSE)


# ======================================
# Création du fichier BED pour bedtools
# ======================================

# Format BED : chr | start (0-based) | end | ID (sonde)
creer_bed <- function(dt, fenetre_kb) {
  fenetre <- fenetre_kb * 1000

  bed <- data.table(
    chr = dt$chr,
    start = pmax(0, dt$start_position - fenetre),
    end = dt$end_position + fenetre,
    ID = paste0(dt$probe_id, "__", dt$symbol)
  )

  # Trier les coordonnées
  bed[, chr_num := suppressWarnings(as.numeric(gsub("chr", "", chr)))]
  bed[chr == "chrX", chr_num := 23]
  bed[chr == "chrY", chr_num := 24]
  setorder(bed, chr_num, start)
  bed[, chr_num := NULL]

  return(bed)
}

bed_10kb <- creer_bed(probe_gene, fenetre_kb = 10)
bed_50kb <- creer_bed(probe_gene, fenetre_kb = 50)

fwrite(bed_10kb, "results/region_cis/regions_cis_10kb.bed", sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(bed_50kb, "results/region_cis/regions_cis_50kb.bed", sep = "\t", col.names = FALSE, quote = FALSE)