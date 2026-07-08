# Projet de stage -- Expression de variants rares

## Objectifs 
- Adapter la méthodologie établie par Mazarotto et al. (2025)
- Détecter les profils d'expression aberrante
- Identifier les variants rares asssociés à une expression génique aberrante
- Quantifier le nombre de gènes affectés par la présence de variants rares
- Tester l'association entre le score IOGC et le statut

## 0. Filtrage des sondes
- Sonde déclarée exprimée si la p-value Illumina Detection Test < 0.05
- Sondes suffisamment exprimées :
        - Dans 75 % et + des individus
        - Dans au moins une famille ayant ≥ 4 membres 

## 1. Détection de l'expression aberrante
### 1.1 Construction du modèle linéaire mixte par sonde
        log_2⁡𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛~ 𝑎𝑔𝑒+𝑠𝑒𝑥𝑒+𝑃𝐶1+ … +𝑃𝐶10 ℎ𝑜𝑢𝑠𝑒𝑘𝑒𝑒𝑝𝑖𝑛𝑔 +1 | 𝑚𝑎𝑡𝑟𝑖𝑐𝑒 𝑑𝑒 𝑘𝑖𝑛𝑠ℎ𝑖𝑝
- PC1...PC10 = capture la variation technique
- Matrice de kinship = correction pour la structure familiale

        log_2⁡𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛~ 𝑎𝑔𝑒+𝑠𝑒𝑥𝑒+𝑃𝐶1+ … +𝑃𝐶10 ℎ𝑜𝑢𝑠𝑒𝑘𝑒𝑒𝑝𝑖𝑛𝑔 +1 | 𝑚𝑎𝑡𝑟𝑖𝑐𝑒 𝑑𝑒 𝑘𝑖𝑛𝑠ℎ𝑖𝑝
- Ajout du modèle linéaire mixte -> évaluer l'impact de l'effet du lead cis-eQTL

* lead cis-eQTL 
            𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 𝑎𝑗𝑢𝑠𝑡é𝑒 ~ 𝑔é𝑛𝑜𝑡𝑦𝑝𝑒 (0 | 1 | 2)+𝑝ℎé𝑛𝑜𝑡𝑦𝑝𝑒 (𝑆𝑍 ┤|𝐵𝐷 | 𝐶𝐿) 
### 1.2 Calcul des scores Z par sondes
        𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 𝑎𝑗𝑢𝑠𝑡é𝑒(𝑠𝑜𝑛𝑑𝑒 𝑖, 𝑖𝑛𝑑𝑖𝑣𝑖𝑑𝑢 𝑗)=  log_2⁡〖𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛(𝑖, 𝑗)  −𝑣𝑎𝑙𝑒𝑢𝑟 𝑝𝑟é𝑑𝑖𝑡𝑒(𝑖, 𝑗)〗
- Définition de l'expression "normale" selon les non-atteints
        𝑍(𝑖, 𝑗)=  ((𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 𝑎𝑗𝑢𝑠𝑡é𝑒(𝑖, 𝑗)  −𝑚𝑜𝑦_𝑁𝐴𝑅(𝑖)))/(é𝑐𝑎𝑟𝑡_𝑡𝑦𝑝𝑒_𝑁𝐴𝑅(𝑖))
- Définition de l'expression "normale" selon les individus atteints (seulement appliqué chez les individus atteints)
        𝑍(𝑖, 𝑗)=  ((𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 𝑎𝑗𝑢𝑠𝑡é𝑒(𝑖, 𝑗)  −𝑚𝑜𝑦_𝐴(𝑖)))/(é𝑐𝑎𝑟𝑡_𝑡𝑦𝑝𝑒_𝐴(𝑖))
* Évaluer si le groupe de référence influence le nombre de outliers identifiés chez les atteints
### 1.3 Identification de l'expression aberrante
Les sondes dépassant le seuil d'expression aberrante **|Z|≥2** sont considérées comme étant aberrante
Les individus avec un nombre de sondes dépassant le seuil adaptif **Q3 + 1.5 * IQR** sont retirés

## 2. Identification des variants rares associés aux outliers
### 2.1 Définition des régions autour des gènes
Obtention des coordonnées des sondes
    - Pour chaque sonde, définir les bornes du gènes correspondant
Définition de la fenêtre
    - ±10 kb
    - ±50 kb
Génération d'un fichier avec les variants rares "candidats"
    - Intersection avec les variants rares
### 2.2 Filtrage des variants dits "candidats" associés aux outliers
Traitement des données
    - Nettoyage et exclusion des symboles manquants
Critère de présence
    - Retenir uniquement les variants portés par au moins un individu identifié comme outlier
Critère de direction cohérente
    - Tous les individus ayant le même variant doivent présenter l'expression dans le même sens

## 3. Calcul du burden de variants rares
### 3.1 Calcul du score IOGC par individu
Score IOGC = nombre de gènes distincts pour lesquelles un individu porte au moins un variant
    - Plusieurs variants dans le même gène compte pour 1
    - Le score n'est pas lié au phénotype
### 3.2 Régression logistique selon le statut
        𝑝ℎé𝑛𝑜𝑡𝑦𝑝𝑒 ~ 𝐼𝑂𝐺𝐶+â𝑔𝑒+𝑠𝑒𝑥𝑒+𝑃𝐶1+ … + 𝑃𝐶𝑘

## 4. Vérification de la robustesse des résultats

## Références 

Mazzarotto F, Gennarelli M, Murray G, Osimo EF. The complementary roles of rare variant burden scores and common variant polygenic risk scores in genetic risk prediction of complex disorders. bioRxiv, 2025. doi: 10.1101/2025.09.15.676308

Chagnon YC, Maziade M, Paccalet T, Croteau J, Fournier A, Roy MA, Bureau A. A multimodal attempt to follow-up linkage regions using RNA expression, SNPs and CpG methylation in schizophrenia and bipolar disorder kindreds. European Journal of Human Genetics, 28:499–507, 2020. doi: 10.1038/s41431-019-0526-y

Ferraro NM, Strober BJ, Einson J, et al. Transcriptomic signatures across human tissues identify functional rare genetic variation. Science, 369(6509):eaaz5900, 2020. doi: 10.1126/science.aaz5900

Smail C, Ferraro NM, Hui Q, et al. Integration of rare expression outlier-associated variants improves polygenic risk prediction. The American Journal of Human Genetics, 109:1055–1064, 2022. doi: 10.1016/j.ajhg.2022.04.015

Ongen H, Buil A, Brown AA, Dermitzakis ET, Delaneau O. Fast and efficient QTL mapper for thousands of molecular phenotypes. Bioinformatics, 32(10):1479–1485, 2016. doi: 10.1093/bioinformatics/btv722