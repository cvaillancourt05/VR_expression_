# Étapes du projet


## 1. Détection de l'expression aberrante

### Prétraitement des données (expression ajustée)
Intensités des sondes déjà corrigées et normalisées
Filtrage des sondes
  - Sonde déclarée exprimée si p-value Illumina detection test < 0.05
  - Sondes suffisamment exprimées : retenu si exprimées dans 75+% et dans au moins une famille >= 4 membres


### 1.1 Calcul des scores Z par sonde

#### a) Construction du modèle linéaire mixte par sonde

    log2_expression ~ age + sexe + PC1..PC10 housekeeping - (lead cis-eQTN) + (1/ matrice kinship)

 - PC1 à PC10 (housekeeping): capture la variation technique
 - Matrice kinship : correction pour la stucture familiale des données
 - lead cis-eQTN : isoler la variation due aux variants rares
   
    → Comparaison avec/sans retrait du lead eQTN pour évaluer l'impact
   
    → Pour chaque sonde, balayage des SNPs dans +/- 1 Mb, sélection de celui avec la plus petite p-value (lead) (Smail)
   
	Selon GTEx:
        cis-eQTL : SNP situé dans une fenêtre de ±1 Mb autour du site d’initiation de la transcription (TSS) du gène
	      lead cis-eQTL :  présente la plus petite p‑value d’association avec l’expression du gène


#### b) Extraction des résidus

    résidu(sonde i, individu j) = log2_expression(i,j) - valeur_prédite(i,j)


#### c) Centrage et réduction (moyenne et écart-type)

Pour définir l'expression "normale" selon les individus non-atteints:

    Z(i,j) = ( résidu(i,j)- moy_nar(i) ) / SD_nar(i)

Pour définir l'expression "normale" selon les individus atteints (seulement applicable chez les atteints) :

    Z(i,j) = ( résidu(i,j)- moy_cas(i) ) / SD_cas(i)

→ Comparaison des deux approches pour évaluer si le groupe de référence influence le nombre de outliers identifiés chez les atteints


### 1.2 Identification des outliers globaux

Seuil d'expression aberrante : |Z|>= 2

  - Retirer les individus avec >1300 sondes dépassant ce seuil (adapté pour le contexte des données)

  - Recalculer les scores Z sur la matrice nettoyée



## 2. Identification des variants rares associés aux outliers (expression)

### 2.1 Vérification de l'assemblage génomique

Effectuer un liftover de hg19 à hg38 (fait dans Ferraro)
  - Pour les sondes d'expression

### 2.2 Définir les régions autour des gènes

#### a) Obtenir les coordonnées des sondes

Pour chaque sonde, définir les bornes du gène correspondant au gène

Avec biomaRt à partir des noms des gènes

#### b) Définir la fenêtre
Tester les deux fenêtres et comparer les résultats :

Fenêtre 1 (Mazzarotto) : +/- 10kb autour du gène (plus conservative, plus spécifique)

Fenêtre 2 (Chagnon) : +/- 50kb autour du gène (capture plus de variants, augmente le bruit)



#### c) Création du fichier BED des régions cis

- Générer un fichier BED avec une ligne par gène (coordonnées en hg38)
		
    chr | start_gene - fenetre | end_gene + fenetre | ID

- Intersecter avec les variants rares pour obtenir une liste des variants candidats (extrait identité des variants + gène associé)

  bedtools intersect -> liste des variants dans les régions

- Extraire les génotypes des variants pour tous individus avec plink ou bcftools ( matrice individus x variants avec génotypes codés 0/1/2)

  Avec l'outil PLINK: 
			- Produit un fichier avec génotypes codés en nombre de copies de l'allèle mineur par individu

- Utiliser bedtools pour extraire les variants génotypés dans les régions

(.ped contient deux colonnes identifier fichier .map voir car variant rare faut identifier allèle mineur)
.frq = frequence des allèles -- identiquer frequence pour chaque allèle 



### 2.3 Filtrer les variants "candidats" associés aux outliers

#### a) Critère de présence

- Retenir uniquement les variants portés par au moins un individu identifié outlier (|Z|>= 2)

- Un variant présent chez un non-outlier pour ce même gène est exclu

#### b) Critère de direction cohérente

Si un variant est porté par plusieurs individus outliers:

- tous doivent présenter l'expression aberrante dans le même sens

- si ce sont des directions opposées, le variant est exclu


#### Résultat : liste de variants candidats avec

  - identifiant du variant

  - sonde associée

  - individus outlier

  - Direction de l'expression aberrante

  - MAF


## 3. Calcul du burden de variants rares (IOGC)

### 3.1 Calcul du score IOGC par individu

IOGC = nombre de gènes distincts pour lesquels un individu porte au moins un variant candidat

  - Plusieurs variants dans le même gène comptent pour 1

  - Le score n'est pas lié au phénotype


Version simplifiée de Mazzarotto -- comptage brut sans pondération par la direction

### 3.2 Association du score IOGC avec les phénotypes (SZ, BD, CL)

Modèle de régression logistique à effets mixtes pour chaque phénotype :

      phénotype ~ IOGC + âge + sexe + PC1..PCk_génotype + (1 / kinship)

- phénotype : variable binaire, donc testé séparément
    SZ vs NAR; BD vs NAR; CL vs NAR

- PC1..PCk_génotype : composantes principales de la matrice de génotype

- (1 / kinship) : pour respecter la dépendance des sujets à cause des familles


Produit un coefficient

- Positif significatif = plus haute probabilité d'être cas en fonction des prédicteurs

Implémenté avec GMMAT

### 3.3 Vérification de la robustesse

#### a) Permutation des scores IOGC

Échanger aléatoirement les scores IOGC entre familles

- Brise le lien entre le score et le phénotype en préservant la structure intra-famille

Répéter N fois

- Si les associations disparaissent avec les scores permutés -> les résultats sont robustes

#### b) Analyse de sensibilité sur les paramètres

Faire varier le seuil de Z

  - Seuils plus stricts donnent des variants probablement plus causaux

Faire varier la fenêtre génomique

  - évalue si les associations sont stables