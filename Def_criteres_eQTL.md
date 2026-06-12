# Définition et critères des cis-eQTL selon GTEx

cis-eQTL : SNP situé dans une fenêtre de ±1 Mb autour du site d’initiation de la transcription (TSS) du gène

cis-eQTL significatif : p‑value après correction pour tests multiples (par permutation ou Bonferroni adaptatif, contrôlant le false discovery rate à 5 %) soit inférieure au seuil calculé pour ce gène

lead cis-eQTL :  présente la plus petite p‑value d’association avec l’expression du gène

# Méthode GTEx pour identifier les cis-eQTL

GTEx utilise un outil pour faire le mapping des cis-eQTL --> FastQTL
		- FastQTL n'est plus maintenu. TensorQTL est le nouveau outil, mais il est fonctionnellement équivalent. 

### Fenêtre cis : +/- 1Mb autour du site de début de transcription (TSS)

### Test d'association - régression linéaire

Pour chaque paire gène-variant :

		expression_résiduelle ~ génotype (0/1/2) + phénotype (SZ/BD/CL) + covariables


### Correction pour tests multiples

Correction intra-gène (avec TensorQTL)

* Permutations de l'expression (environ une centaine)
* Récupération de la meilleure p-value par permutation
* Modélisation de la distribution par une loi beta pour obtenir des p-value ajustées

--Potentiellement computationnellement lourd -> GPU

Correction genome-wide --correction FDR

* Seuil de <5% sur les p-value pour être considéré cis-eQTL significatif

#### d) Lead cis-eQTL par gène

Variant avec la p-value d'association la plus faible dans la fenêtre

# Pipeline d'analyse cis-eQTL - implémentation

### Préparation des entrées via prep_tensorQTL.R

1. Fusion de la matrice d'expression avec les coordonnées en hg38
2. Filtrage - individus au statut inconnu retirés
3. Formatage BED - 
		#Chr, start, end, ID
4. Construction du fichier phénotypes
5. Compression bgzf + indexation Tabix

Sorties : 
* Matrice d'expression au format BED
* Covariable diagnostique

### Sorties de TensorQTL

* Les p-values nominales et ajustées pour toutes les paires gène-variant (tensorqtl.cis_qtl.txt)