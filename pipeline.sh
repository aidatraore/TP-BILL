####################### Pour la qualité du sequençage, on utilise pycoQc est un outil utilisé pour analyser et visualiser la qualité des données issues de sequençage oxford Nanapore. 
srun pycoQC -f  /students/BILL/2026-BILL/P25_P27Aprem/sequencing_summary.txt -o sequencing.html

#### on obtient un fichier.html  avec le pycoQC pour visualiser les resultats de la qualité du séquençage############

############### pipeline d'analyse ##################

######### on va verifier la qualité de nos fichier VCF ###########
bcftools stats ./Inputs/P25/P25-1.trimed1000.sv_sniffles.vcf > P25-1_stats.txt

#### visualiser les stats

less P25-1_stats.txt

########## Résumé graphique automatique

plot-vcfstats P25-1_stats.txt  -p stats_plots


######### pour le filtrage

### faire le filtre directement Filtrer sur DP depuis INFO avec bcftools

bcftools filter -i 'QUAL>=30 && INFO/DP>=10' mon_fichier.vcf -o mon_fichier_filtered.vcf

####### vu que on n'a pas de DP donc le filtrage a été fait sans le DP

######## une fois le filtrage on merge les fivhier vcf avec bcftools merge.

#### pour merger, on foit compresser et indexer tout les fichier vcf

# Compresser
bgzip P25_filtered.vcf
bgzip P27_filtered.vcf

# Indexer
bcftools index P25_filtered.vcf.gz
bcftools index P27_filtered.vcf.gz

## fusionner les VCF
bcftools merge P25_filtered.vcf.gz P27_filtered.vcf.gz -o merged.vcf.gz -O z


### fusion en intersection
bcftools isec -n=2 P25_filtered.vcf.gz P27_filtered.vcf.gz -w1 -O z -o common.vcf.gz


### extraction des infos dans le vcf


for f in P25/*.snp.vcf
do
  sample=$(basename "$f" .snp.vcf)
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO\t%FORMAT\n' "$f" > "Analyses/P25/${sample}_table_snp.vcf"
done


#### ### pipeline

Filtrage QUAL
      ↓
Merge sous-cultures par passage
      ↓
Rename samples
      ↓
Merge global temporel
      ↓
Extraction table
      ↓
Analyse R / visualisation








