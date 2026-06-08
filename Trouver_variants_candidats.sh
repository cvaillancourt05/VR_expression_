#!/bin/bash

#SBATCH --account=def-bureau
#SBATCH --mem=8G
#SBATCH --time=02:00:00

module load plink/2.00-20231024-avx2
module load bcftools/1.22
module load bedtools/2.31.0

# ------------------------------------------------------------------------------------------------------
# Trouver_variants_candidats.sh
# Identifie les variants rares dans les régions cis des sondes d'expression et extrait leurs génotypes
# 
# Entrées :
#  - VCF compressés par chromosome
#  - Régions cis par fenêtre -- region_cis_10kb/50kb.bed
#
# Sorties : 
#  - Paires variant-sonde -- variants_dans_regions_<fenetre>.bed
#  - Liste des variants dans les régions cis -- liste_variants_<fenetre>.txt
#  - Fréquences alléliques des variants retenus (pour annotation) -- freq_variants_<fenetre>.frq 
# ------------------------------------------------------------------------------------------------------

# --------
# Chemins
# --------

# --Répertoire des VCF par chromosome
seq_dir="/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
# --Préfixe commun des VCF
vcf_prefix="impute5_gigi2_combined_seq_RV"

# --Préfixe commun des fichiers
frq_prefix="impute5_gigi2_combined_seq_RV_FINAL"
# --Répertoire des fichiers frq
frq_dir="/lustre09/project/6033529/schizo/data/WGS_bs_2022/freq_RV/"

# --Répertoire de sortie pour les fichiers intermédiaires par chromosome
out_dir="/home/chloev/links/projects/def-bureau/chloev/liste_variants/sorties"
# --Répertoire  contenant les fichiers de régions cis
misc_dir="/home/chloev/links/projects/def-bureau/chloev/liste_variants"
# --Répertoire de sortie pour les fichiers fusionnés
merge_dir="${out_dir}/merge"

# --------------------------
# Traitement par chromosome
# --------------------------

for chr in $(seq 1 22); do

    vcf="${seq_dir}/${vcf_prefix}_chr${chr}.vcf.gz"

    # --Extraction des coordonnées du VCF au format BED (CHROM, POS0, END, ID)
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$vcf" \
    | sed 's/^chr//' \
    | sort -k1,1 -k2,2n > "${out_dir}/variants_chr${chr}.bed"

    for FENETRE in 10kb 50kb; do
        
        # --Filtrage de la région cis pour le chromosome actuel
        grep -P "^chr${chr}\t" "${misc_dir}/regions_cis_${FENETRE}.bed" \
        | sed 's/^chr//' > "${out_dir}/regions_${FENETRE}_chr${chr}.bed"

        # --Intersection variants x régions cis
        bedtools intersect -a "${out_dir}/variants_chr${chr}.bed" -b "${out_dir}/regions_${FENETRE}_chr${chr}.bed" -wa -wb | \
        awk 'BEGIN{OFS="\t"} { print $4, $8 }' | sort -u > "${out_dir}/variants_dans_regions_${FENETRE}_chr${chr}.bed"

        # --Liste des variants à extraire
        cut -f1 "${out_dir}/variants_dans_regions_${FENETRE}_chr${chr}.bed" | sort -u  > "${out_dir}/liste_variants_${FENETRE}_chr${chr}.txt"

    done
done
     
# ---------------------------
# Fusion globale par fenêtre
# ---------------------------

for FENETRE in 10kb 50kb; do

    # --Fusion des paires variant-sonde de tous les chromosomes
    cat "${out_dir}"/variants_dans_regions_${FENETRE}_chr*.bed | sort -u > "${merge_dir}/variants_dans_regions_${FENETRE}.bed"

    # --Fusion des listes de variants uniques
    cat "${out_dir}"/liste_variants_${FENETRE}_chr*.txt | sort -u > "${merge_dir}/liste_variants_${FENETRE}.txt"

    # --Extraction des fréquences alléliques pour annotation seulement
    head -1 "${frq_dir}/${frq_prefix}_chr1.frq.cc" > "${merge_dir}/freq_variants_${FENETRE}.frq"

    for chr in $(seq 1 22); do
        frq="${frq_dir}/${frq_prefix}_chr${chr}.frq.cc"
        liste="${merge_dir}/liste_variants_${FENETRE}.txt"
        awk 'NR==FNR {ids[$1]=1; next} FNR > 1 && $2 in ids' "$liste" "$frq" >> "${merge_dir}/freq_variants_${FENETRE}.frq"
    done 

    # --Extraction finale des génotypes par chromosome depuis les VCF
    for chr in $(seq 1 22); do
        bcftools view --include ID=@"${merge_dir}/liste_variants_${FENETRE}.txt" "${seq_dir}/${vcf_prefix}_chr${chr}.vcf.gz" \
        | bcftools query --include 'GT!="0/0" & GT!="0|0" & GT!="./."' -f '[%SAMPLE\t%ID\t%GT\n]' \
        | awk 'BEGIN{OFS="\t"} {
            gt = $3
            if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") geno=1
            else if (gt=="1/1" || gt=="1|1") geno=2
            else next
            print $1, $2, geno
          }' \
        > "${merge_dir}/porteurs_${FENETRE}_chr${chr}.tsv"
    done    
done


rm "${out_dir}"/*_chr*.bed "${out_dir}"/*_chr*.txt "${out_dir}"/*_chr*.log