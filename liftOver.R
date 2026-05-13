#Liftover -- hg19 à hg38

library(GenomicRanges)
library(rtracklayer)

# Lecture du fichier
annotation <- read.delim("data/Illumina HumanHT-12 V4.0 expression beadchip/HumanHT-12_V4_0_R2_15002873_B.txt", skip = 8, stringsAsFactors = FALSE)

# Récupérer les sondes humaines
annotation <- annotation[annotation$Species == "Homo sapiens", ]
annotation <- annotation[as.character(annotation$Chromosome) %in% c(as.character(1:22), "X", "Y"), ]

# Parsing des coordonnées 
parse_probe <- function(chr, coord_string, probe_id, illumina_id, symbol) {

    pieces <- unlist(strsplit(coord_string, ":"))
    ranges <- lapply(pieces, function(x) {

        pos <- unlist(strsplit(x, "-"))
        data.frame(
            chr = paste0("chr", chr),
            start = as.numeric(pos[1]),
            end = as.numeric(pos[2]),
            probe_id = probe_id,
            illumina_id = illumina_id,
            symbol = symbol,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, ranges)
}

# Parsing de toutes les sondes
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

# Retirer les chromosomes invalides
all_ranges <- all_ranges[
    gsub("chr", "", all_ranges$chr) %in%
        c(as.character(1:22), "X", "Y"),
]

# Construction des GRanges hg19
gr_hg19 <- GRanges(seqnames = all_ranges$chr, ranges = IRanges(start = all_ranges$start, end = all_ranges$end), probe_id = all_ranges$probe_id, illumina_id = all_ranges$illumina_id, symbol = all_ranges$symbol)

# Téléchargement de la chain file
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
chain <- import.chain("data/Illumina HumanHT-12 V4.0 expression beadchip/hg19ToHg38.over.chain")

# LiftOver
gr_hg38_list <- liftOver(gr_hg19, chain)
keep <- elementNROWS(gr_hg38_list) == 1
gr_hg38 <- unlist(gr_hg38_list[keep])

# Construction de la table finale
lifted <- data.frame(
  probe_id = gr_hg38$probe_id,
  illumina_id = gr_hg38$illumina_id,
  symbol = gr_hg38$symbol,
  chr_hg38 = as.character(seqnames(gr_hg38)),
  start_hg38 = start(gr_hg38),
  end_hg38 = end(gr_hg38),
  stringsAsFactors = FALSE
)

# Sauvegarde
#write.csv(lifted, "data/Illumina HumanHT-12 V4.0 expression beadchip/resultat/HumanHT12v4_hg38_annotations.csv", row.names = FALSE)