####################### Pour la qualité du sequençage, on utilise pycoQc est un outil utilisé pour analyser et visualiser la qualité des données issues de sequençage oxford Nanapore. 
srun pycoQC -f  /students/BILL/2026-BILL/P25_P27Aprem/sequencing_summary.txt -o sequencing.html

#### on obtient un fichier.html  avec le pycoQC pour visualiser les resultats de la qualité du séquençage############

############### pipeline d'analyse ##################

######### on va verifier la qualité de nos fichier VCF ###########
bcftools stats ./Inputs/P25/P25-1.trimed1000.sv_sniffles.vcf > P25-1_stats.txt
#### stat pour touts les replicats
for f in P25/*.sv_sniffles.vcf
do
bcftools stats "$f" > "Analyses/$(basename ${f%.vcf})_stats_sv_snif.txt"
done

for f in P27/*.snp.vcf
do
bcftools stats "$f" > "Analyses/$(basename ${f%.vcf})_stats_snp.txt"
don



#### visualiser les stats

less P25-1_stats.txt

########## Résumé graphique automatique

plot-vcfstats P25-1_stats.txt  -p stats_plots


######### pour le filtrage

### faire le filtre directement Filtrer sur DP depuis INFO avec bcftools

bcftools filter -i 'QUAL>=30 && INFO/DP>=10' mon_fichier.vcf -o mon_fichier_filtered.vcf

####### vu que on n'a pas de DP donc le filtrage a été fait sans le DP
###### stat pour tout les vcf filtré
for f in *snp-filtr
do
bcftools stats "$f" > "$(basename ${f%.vcf})_stats_snp.txt"
done

######## une fois le filtrage on merge les fichiers vcf avec bcftools merge.

#### pour merger, il faut compresser et indexer tout les fichier vcf

###### Compresser
bgzip P25_filtered.vcf
bgzip P27_filtered.vcf

# Indexer
bcftools index P25_filtered.vcf.gz
bcftools index P27_filtered.vcf.gz

## fusionner les VCF
bcftools merge P25_filtered.vcf.gz P27_filtered.vcf.gz -o merged.vcf.gz -O z


### fusion en intersection
bcftools isec -n=2 P25_filtered.vcf.gz P27_filtered.vcf.gz -w1 -O z -o common.vcf.gz

#### renommer toute les variables à partir de de chrom de ba base mergé
## renommer

sed -i 's/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\t2:SAMPLE\t3:SAMPLE\t4:SAMPLE\t5:SAMPLE\t6:SAMPLE/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP15_2\tP25_2\tP27_2\tP30_2\tP50_2\tP90_2/' basemergefr.vcf


### extraction des infos dans le vcf


for f in P25/*.snp.vcf
do
  sample=$(basename "$f" .snp.vcf)
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO\t%FORMAT\n' "$f" > "Analyses/P25/${sample}_table_snp.vcf"
done

#### calculer le depth a partir de samtools
for f in P27/*aligned.sorted.bam
do
  sample=$(basename "$f" .aligned.sorted.bam)
  samtools depth -a -q 20 -Q 20 -F 256 -F 1024 "$f" > "Analyses/P27/${sample}_depth.txt"
  echo "Depth pour $sample généré."
done

for f in P25/*aligned.sorted.bam
do
  sample=$(basename "$f" .aligned.sorted.bam)
  samtools depth -a -q 20 -Q 20 -F 256 -F 1024 "$f" > "Analyses/P25/${sample}_depth.txt"
  echo "Depth pour $sample généré."
done

###### extraire les infos dans le fichier flagstats
echo -e "Sample\tTotal_reads\tPrimary_reads\tSecondary\tSupplementary\tDuplicates\tMapped\tPrimary_mapped" > summary_flagstat_P25.tsv

output="Analyses/P25/summary_flagstat_P25.txt"

# créer entête
echo -e "Sample\tTotal\tPrimary\tSecondary\tSupplementary\tDuplicates\tMapped\tPrimary_mapped" > "$output"

for f in P25/*.flagstat
do
    sample=$(basename "$f" .flagstat)

    total=$(grep "in total" "$f" | awk '{print $1}')
    primary=$(grep "primary$" "$f" | awk '{print $1}')
    secondary=$(grep "secondary" "$f" | awk '{print $1}')
    supplementary=$(grep "supplementary" "$f" | awk '{print $1}')
    duplicates=$(grep "duplicates" "$f" | head -1 | awk '{print $1}')
    mapped=$(grep " mapped (" "$f" | head -1 | awk '{print $1}')
    rate=$(grep " mapped (" "$f" | head -1 | awk '{print $2}')

    echo -e "$sample\t$total\t$primary\t$secondary\t$supplementary\t$duplicates\t$mapped\t$primary_mapped" >> "$output"

done

output="Analyses/P27/summary_flagstat_P27.txt"
# créer entête
echo -e "Sample\tTotal\tPrimary\tSecondary\tSupplementary\tDuplicates\tMapped\tPrimary_mapped" > "$output"

for f in P27/*.flagstat
do
    sample=$(basename "$f" .flagstat)

    total=$(grep "in total" "$f" | awk '{print $1}')
    primary=$(grep "primary$" "$f" | awk '{print $1}')
    secondary=$(grep "secondary" "$f" | awk '{print $1}')
    supplementary=$(grep "supplementary" "$f" | awk '{print $1}')
    duplicates=$(grep "duplicates" "$f" | head -1 | awk '{print $1}')
    mapped=$(grep " mapped (" "$f" | head -1 | awk '{print $1}')
    rate=$(grep " mapped (" "$f" | head -1 | awk '{print $2}')

    echo -e "$sample\t$total\t$primary\t$secondary\t$supplementary\t$duplicates\t$mapped\t$primary_mapped" >> "$output"

done


bcftools merge $(ls P15-[1-5]*trimed1000.snp_filtered.vcf | grep -v barcode) -o P15_1-5_merged.vcf -O v


for f in *trimed1000.snp_filtered.vcf; do
    # ignorer les fichiers contenant 'barcode'
    if [[ "$f" != *barcode* ]]; then
        echo "Compression et indexation de $f..."
        bgzip -c "$f" > "${f}.gz"
        tabix -p vcf "${f}.gz"
    fi
done


bcftools merge P15-[1-5]*trimed1000.snp_filtered.vcf.gz -O v -o P15_1-5_merged.vcf
echo "Fusion P15-1 à P15-5 terminée : P15_1-5_merged.vcf"


bcftools merge P15-[1-5]*trimed1000.snp_filtered.vcf.gz -O v -o P15_1-5_merged.vcf

bcftools merge P15


--force-samples


bcftools merge --force-samples  P90-1_merged.trimed1000.snp_filtered.vcf.gz P90-2_merged.trimed1000.snp_filtered.vcf.gz P90-3_merged.trimed1000.snp_filtered.vcf.gz P90-4_merged.trimed1000.snp_filtered.vcf.gz  P90-5_merged.trimed1000.snp_filtered.vcf.gz -O v -o P90_1-5_merged.vcf


bcftools merge --force-samples  P90-6_merged.trimed1000.snp_filtered.vcf.gz P90-7_merged.trimed1000.snp_filtered.vcf.gz P90-8_merged.trimed1000.snp_filtered.vcf.gz P90-9_merged.trimed1000.snp_filtered.vcf.gz P90-10_merged.trimed1000.snp_filtered.vcf.gz -O v -o P90_6-10_merged.vcf



## renommé
echo -e "SAMPLE\tP90-1" > mapping.txt
echo -e "2:SAMPLE\tP90-2" >> mapping.txt
echo -e "3:SAMPLE\tP90-3" >> mapping.txt
echo -e "4:SAMPLE\tP90-4" >> mapping.txt
echo -e "5:SAMPLE\tP90-5" >> mapping.txt


bcftools reheader -s mapping.txt -o P90_1-5_renamed.vcf P90_1-5_merged.vcf

## verification

bcftools query -l P90_1-5_renamed.vcf


## renommé
echo -e "SAMPLE\tP90-6" > mapping.txt
echo -e "2:SAMPLE\tP90-7" >> mapping.txt
echo -e "3:SAMPLE\tP90-8" >> mapping.txt
echo -e "4:SAMPLE\tP90-9" >> mapping.txt
echo -e "5:SAMPLE\tP90-10" >> mapping.txt


bcftools reheader -s mapping.txt -o P90_6-10_renamed.vcf P90_6-10_merged.vcf

## verification

bcftools query -l P90_6-10_renamed.vcf

for f in *_renamed.vcf
do
    bgzip -c $f > ${f}.gz
    tabix -p vcf ${f}.gz
done

## merge
bcftools merge P*_merged.vcf.gz -O v -o ALL_PASSAGES.vcf

bcftools merge P15_1-5_renamed.vcf.gz P25_1-5_renamed.vcf.gz  P27_1-5_renamed.vcf.gz  P30_1-5_renamed.vcf.gz P50_1-5_renamed.vcf.gz P65_1-5_renamed.vcf.gz P90_1-5_renamed.vcf.gz -O v -o ALL_PASSAGES.vcf

bcftools merge P15_1-5_renamed.vcf.gz P25_1-5_renamed.vcf.gz  P27_1-5_renamed.vcf.gz  P30_1-5_renamed.vcf.gz P50_1-5_renamed.vcf.gz P65_1-5_renamed.vcf.gz P90_1-5_renamed.vcf.gz -O v -o Base_froid.vcf

bcftools merge P15_6-10_renamed.vcf.gz P25_6-10_renamed.vcf.gz P27_6-10_renamed.vcf.gz P30_6-10_renamed.vcf.gz P50_6-10_renamed.vcf.gz P65_6-10_renamed.vcf.gz
P90_6-10_renamed.vcf.gz -O v -o Base_chaud.vcf

#### supprimer dans le fichier vcf toutes les lignes commencant par ## afin d'avoir une base .tsv froi et chaud
grep -v '^##' Base_froid.vcf | sed 's/^#//g' | tr '\t' '\t' > Base_froid.tsv

grep -v '^##' Base_chaud.vcf | sed 's/^#//g' | tr '\t' '\t' > Base_chaud.tsv








