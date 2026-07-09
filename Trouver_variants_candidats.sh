#!/bin/bash

#SBATCH --account=def-bureau
#SBATCH --mem=8G
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4

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

# --Répertoires + préfix
seq_dir="/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
vcf_prefix="impute5_gigi2_combined_seq_RV"
frq_prefix="impute5_gigi2_combined_seq_RV_FINAL"
frq_dir="/lustre09/project/6033529/schizo/data/WGS_bs_2022/freq_RV/"
out_dir="/home/chloev/links/projects/def-bureau/chloev/liste_variants/sorties"
misc_dir="/home/chloev/links/projects/def-bureau/chloev/liste_variants"
merge_dir="${out_dir}/merge"

# --------------------------
# Traitement par chromosome
# --------------------------

for chr in $(seq 1 22); do

    vcf="${seq_dir}/${vcf_prefix}_chr${chr}.vcf.gz"

    # --Extraction des coordonnées du VCF (CHROM, POS, ID, AF)
    # ---si AF n'est pas dans l'INFO, on utilise le plugin fill-tags
    bcftools +fill-tags "$vcf" -- -t AF | bcftools query -f '%CHROM\t%POS0\t%END\t%ID\t%INFO/AF\n' \
    | sed 's/^chr//' > "${out_dir}/variants_chr${chr}.txt"

    for FENETRE in 10kb 50kb; do
        
        # --Filtrage de la région cis pour le chromosome actuel
        grep -P "^chr${chr}\t" "${misc_dir}/regions_cis_${FENETRE}.bed" \
        | sed 's/^chr//' > "${out_dir}/regions_${FENETRE}_chr${chr}.bed"

        # --Intersection variants x régions cis
        bedtools intersect -a "${out_dir}/variants_chr${chr}.txt" -b "${out_dir}/regions_${FENETRE}_chr${chr}.bed" -wa -wb \
        | awk 'BEGIN{OFS="\t"} { print $4, $9, $5 }' > "${out_dir}/variants_dans_regions_${FENETRE}_chr${chr}.txt"

    done
done
     
# ---------------------------
# Fusion globale par fenêtre
# ---------------------------

for FENETRE in 10kb 50kb; do

    # --Fusion des paires variant-sonde de tous les chromosomes
    cat "${out_dir}"/variants_dans_regions_${FENETRE}_chr*.txt > "${merge_dir}/temp_fusion_${FENETRE}.txt"
    # Créer la liste des paires variant-sonde (ton format .bed original)
    awk 'BEGIN{OFS="\t"} {print $1, $2}' "${merge_dir}/temp_fusion_${FENETRE}.txt" > "${merge_dir}/variants_dans_regions_${FENETRE}.bed"
    # Liste unique des variants
    awk '{print $1}' "${merge_dir}/temp_fusion_${FENETRE}.txt" | sort -u > "${merge_dir}/liste_variants_${FENETRE}.txt"
    # Liste des variants à inverser (AF > 0.5)
    awk '$3 > 0.5 {print $1}' "${merge_dir}/temp_fusion_${FENETRE}.txt" | sort -u > "${merge_dir}/variants_inverser_${FENETRE}.txt"

    # --Extraction finale des génotypes par chromosome depuis les VCF
    for chr in $(seq 1 22); do
        bcftools view --threads 4 --include ID=@"${merge_dir}/liste_variants_${FENETRE}.txt" "${seq_dir}/${vcf_prefix}_chr${chr}.vcf.gz" \
        | bcftools query --include 'GT!="./."' -f '[%SAMPLE\t%ID\t%GT\n]' \
        | awk 'BEGIN{OFS="\t"} {
            if (NR==FNR) { inv[$1]=1; next }
    
                gt = $3
                if (gt=="0/0" || gt=="0|0") alt_count=0
                else if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") alt_count=1
                else if (gt=="1/1" || gt=="1|1") alt_count=2
                else next

                if ($2 in inv) geno = 2 - alt_count
                else geno = alt_count

                if (geno > 0) print $1, $2, geno
                }' "${merge_dir}/variants_inverser_${FENETRE}.txt" - > "${merge_dir}/porteurs_${FENETRE}_chr${chr}.tsv"
    done    
done


rm "${out_dir}"/*_chr*.bed "${out_dir}"/*_chr*.txt "${out_dir}"/*_chr*.log 