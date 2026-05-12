# Définition et critères des cis-eQTL selon GTEx

cis-eQTL : SNP situé dans une fenêtre de ±1 Mb autour du site d’initiation de la transcription (TSS) du gène

cis-eQTL significatif : p‑value après correction pour tests multiples (par permutation ou Bonferroni adaptatif, contrôlant le false discovery rate à 5 %) soit inférieure au seuil calculé pour ce gène

lead cis-eQTL :  présente la plus petite p‑value d’association avec l’expression du gène

## Méthode GTEx pour identifier les cis-eQTL

GTEx utilise un outil pour faire le mapping des cis-eQTL --> FastQTL

#### a) Fenêtre cis : +/- 1Mb autour du site de début de transcription (TSS)

#### b) Test d'association - régression linéaire

Pour chaque paire gène-variant :

		expression résiduelle/ajustée \~ génotype + covariables



#### c) correction pour tests multiples

Correction intra-gène (avec FastQTL)

* Permutations de l'expression (environ une centaine)
* Récupération de la meilleure p-value par permutation
* Modélisation de la distribution par une loi beta pour obtenir des p-value ajustées
--Potentiellement computationnellement lourd

Correction genome-wide --correction FDR

* Seuil de <5% sur les p-value pour être considéré cis-eQTL significatif

#### d) Lead cis-eQTL par gène

Variant avec la p-value d'association la plus faible dans la fenêtre


## Régression linéaire pour identification des lead cis-eQTN

      expression_résiduelle ~ génotype (0/1/2) + phénotype (0/1/2) + covariables

