library(data.table)
library(readxl)

filtre <- read_xlsx("/lustre09/project/6000443/expression_genes/table_expressed_probes_sans_lame_corr_bg_RMA_REML_GCbr_543sujets.xlsx")
sondes_valides <- unique(na.omit(as.character(filtre[[1]])))

# sauvegarde
write.table(data.frame(probe_id = sondes_valides), "/lustre09/project/6000443/expression_genes/liste_sondes_exprimees.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")