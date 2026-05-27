library(data.table)

# =======================
# Chargement des données
# =======================

mat <- fread("data/sondes_exprimees_fam_p05_75percent_549sujets.csv", sep = ";", header = TRUE, data.table = FALSE)
rownames(mat) <- sub("^X", "", mat[[1]])
mat <- mat[, -1]
mat_bool <- as.matrix(mat) == "TRUE" 

# Sondes exprimées doivent faire partie d'une famille ayant au mons 4 membres
sonde_ok <- rowsums(mat_bool, na.rm = TRUE) >= 1
sondes_filtrees <- rownames(mat_bool)[sondes_ok]

expression <- fread("data/Adjusted_expression_values.txt", select = 1)
colnames(expression)[1] <- "probe_id"
sondes_finales(expression$probe_id, sondes_filtrees)

fwrite(data.frame(probe_id = sondes_finales), "results/liste_sondes_exprimees.txt", col.names = FALSE)