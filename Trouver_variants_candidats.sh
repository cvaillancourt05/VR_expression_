#!/bin/bash

# Chemins pour les répertoire
seq_dir="/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
vcf_prefix="impute5_gigi2_combined_seq_RV_FINAL"
out_dir="/home/chloev/links/projects/def-bureau/chloev/liste_variants/sorties"
misc_dir="/home/chloev/links/projects/def-bureau/chloev/liste_variants"
merge_dir="${out_dir}/merge"


# --------------------------
# Traitement par chromosome
# --------------------------

for chr in $(seq 1 22); do

    vcf="${seq_dir}/${vcf_prefix}_chr${chr}.vcf.gz"

    # Extraction des coordonnées du VCF au format BED
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$vcf" \
    | sort -k1,1 -k2,2n > "${out_dir}/variants_chr${chr}.bed"

    for FENETRE in 10kb 50kb; do
        
        # Filtrage de la région cis pour le chromosome actuel
        grep -P "^chr${chr}\t" "${misc_dir}/regions_cis_${FENETRE}.bed" > "${out_dir}/regions_${FENETRE}_chr${chr}.bed"

        # Intersection
        bedtools intersect -a "${out_dir}/variants_chr${chr}.bed" -b "${out_dir}/regions_${FENETRE}_chr${chr}.bed" -wa -wb | \
        awk 'BEGIN{OFS="\t"} { print $4, $8 }' | sort -u > "${out_dir}/variants_dans_regions_${FENETRE}_chr${chr}.bed"

        # Liste des variants à extraire
        cut -f1 "${out_dir}/variants_dans_regions_${FENETRE}_chr${chr}.bed" | sort -u  > "${out_dir}/liste_variants_${FENETRE}_chr${chr}.txt"

        # Extraction des génotypes
        plink2 --vcf "$vcf" \
            --allow-extra-chr \
            --extract "${out_dir}/liste_variants_${FENETRE}_chr${chr}.txt" \
            --export A \
            --out "${out_dir}/genotypes_${FENETRE}_chr${chr}"

    done
done

       
# ---------------
# Fusion globale
# ---------------

for FENETRE in 10kb 50kb; do

    # Fusion variant-sonde
    cat "${out_dir}"/variants_dans_regions_${FENETRE}_chr*.bed | sort -u > "${merge_dir}/variants_dans_regions_${FENETRE}.bed"

    # Fusion des listes de variants uniques
    cat "${out_dir}"/liste_variants_${FENETRE}_chr*.txt | sort -u > "${merge_dir}/liste_variants_${FENETRE}.txt"

    # Fusion des fichiers de fréquences, filtré sur les variants d'intérêt
    head -1 "${seq_dir}/${vcf_prefix}_chr1.frq" > "${merge_dir}/freq_variants_${FENETRE}.frq"

    for chr in $(seq 1 22); do
        frq="${seq_dir}/${vcf_prefix}_chr${chr}.frq"
        liste="${out_dir}/liste_variants_${FENETRE}_chr${chr}.txt"
        awk 'NR==FNR {ids[$1]=1; next} $2 in ids' "$liste" "$frq" >> "${merge_dir}/freq_variants_${FENETRE}.frq"
    done 

    # Extraction finale des génotypes par chromosome
    for chr in $(seq 1 22); do
        plink2 --vcf "${seq_dir}/${vcf_prefix}_chr${chr}.vcf.gz" \
            --allow-extra-chr \
            --extract "${merge_dir}/liste_variants_${FENETRE}.txt" \
            --export A \
            --out "${merge_dir}/genotypes_${FENETRE}_chr${chr}"
    done    
done


rm "${out_dir}"/*_chr*.bed "${out_dir}"/*_chr*.txt "${out_dir}"/*_chr*.log