# ------------------------------------------------------------------------------------------------------------------
# liftOver.R
# Conversion des coordonnées des sondes Illumina HumanHT-12 v4 de hg19 <a hg38
#
# Entrées : 
#  - Fichier d'annotation Illumina 
#  - Fichier chain UCSC hg19Tohg38.over.chain
# Sortie : 
#  - Table 'lifted' (optionnellement sauvegarder en .csv)  --- probe_id | illumina_id | symbol | chr/start/end hg38
# -------------------------------------------------------------------------------------------------------------------

library(GenomicRanges)
library(rtracklayer)

# ----------------------------------------------
# Lecture et nettoyage de l'annotation Illumina
# ----------------------------------------------

# --Lecture du fichier
annotation <- read.delim("data/Illumina HumanHT-12 V4.0 expression beadchip/HumanHT-12_V4_0_R2_15002873_B.txt", skip = 8, stringsAsFactors = FALSE)

# Récupérer uniquement les sondes humaines (1-22, X, Y)
annotation <- annotation[annotation$Species == "Homo sapiens", ]
annotation <- annotation[as.character(annotation$Chromosome) %in% c(as.character(1:22), "X", "Y"), ]

# ----------------------------------------------
# Parsing des coordonnées génomiques des sondes
# ----------------------------------------------

# --Retourne un data frame avec une ligne par coordonnée valide
parse_probe <- function(chr, coord_string, probe_id, illumina_id, symbol) {
    if (is.na(coord_string) || coord_string == "" || coord_string == "0") {
        return(data.frame())
    }
    
    pieces <- unlist(strsplit(coord_string, ";"))
    
    ranges <- lapply(pieces, function(x) {

        if (grepl(":", x)) {
            x <- unlist(strsplit(x, ":"))[2]
        }
        
        pos <- unlist(strsplit(x, "-"))
        start_val <- as.numeric(pos[1])
        end_val <- as.numeric(pos[2])
        
        # ---si le parsing échoue, on retourne NULL
        if (is.na(start_val) || is.na(end_val)) return(NULL)
        
        data.frame(
            chr = paste0("chr", chr),
            start = start_val,
            end = end_val,
            probe_id = probe_id,
            illumina_id = illumina_id,
            symbol = symbol,
            stringsAsFactors = FALSE
        )
    })
    
    # ---filtrer les NULL avant dfusion
    ranges <- ranges[!sapply(ranges, is.null)]
    if (length(ranges) == 0) return(data.frame())
    
    do.call(rbind, ranges)
}

# --Apppliquer le parsing sur toutes les sondes
all_ranges_list <- lapply(seq_len(nrow(annotation)), function(i) {

    parse_probe(
        chr = annotation$Chromosome[i],
        coord_string = annotation$Probe_Coordinates[i],
        probe_id = annotation$Array_Address_Id[i],
        illumina_id = annotation$ILMN_Gene[i],
        symbol = annotation$Symbol[i]
    )
})

all_ranges <- do.call(rbind, all_ranges_list)

# --Retirer les chromosomes invalides
all_ranges <- all_ranges[
    gsub("chr", "", all_ranges$chr) %in%
        c(as.character(1:22), "X", "Y"),
]

# ----------------------
# LiftOver hg19 -> hg38
# ----------------------

# --Construction des GRanges en coordonnées hg19
gr_hg19 <- GRanges(seqnames = all_ranges$chr, ranges = IRanges(start = all_ranges$start, end = all_ranges$end), probe_id = all_ranges$probe_id, illumina_id = all_ranges$illumina_id, symbol = all_ranges$symbol)

# --Fichier à télécharger depuis :
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
chain <- import.chain("data/Illumina HumanHT-12 V4.0 expression beadchip/hg19ToHg38.over.chain")

gr_hg38_list <- liftOver(gr_hg19, chain)
# --Garder uniquement les sondes avec une correspondance unique en hg38
keep <- elementNROWS(gr_hg38_list) == 1
gr_hg38 <- unlist(gr_hg38_list[keep])

# -------------------------------------------------
# Construction de la table finale et déduplication
# -------------------------------------------------

lifted <- data.frame(
  probe_id = gr_hg38$probe_id,
  illumina_id = gr_hg38$illumina_id,
  symbol = gr_hg38$symbol,
  chr_hg38 = as.character(seqnames(gr_hg38)),
  start_hg38 = start(gr_hg38),
  end_hg38 = end(gr_hg38),
  stringsAsFactors = FALSE
)

# --Déduplication par probe_id
# En cas de sondes dupliquées, garder la première occurrence
lifted <- lifted[order(lifted$probe_id, lifted$start_hg38), ]
lifted <- lifted[!duplicated(lifted$probe_id), ]

# --Sauvegarde
# write.csv(lifted, "results/HumanHT12v4_hg38_annotations.csv", row.names = FALSE)