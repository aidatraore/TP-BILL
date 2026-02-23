### Repertoire de travaille
setwd("~/MASTER BIOINFORMATIQUE/TP_Bill/stats_resume")


{
  library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(magrittr)
library(reshape2)
  }

# Lire le TSV
vcf <- read.delim("Base_froid.tsv", stringsAsFactors = FALSE, check.names = FALSE)

# Vérifier les colonnes
head(vcf)

#### DIFFERENCIÉ les snp ,insert, deletion

vcf$typevariant <-ifelse(nchar(vcf$REF) < nchar(vcf$ALT), "Insertion" , 
                         ifelse (nchar(vcf$REF) > nchar (vcf$ALT), "deletion",
                                 "SNP"))

table_variants <- table(vcf$typevariant)
table_variants

print(vcf[vcf$typevariant == "Insertion", ])

##### metre mes sample en colone
vcf_long <- vcf %>%
  pivot_longer(
    cols = -c(CHROM:FORMAT, typevariant),
    names_to = "Sample",
    values_to = "Genotype"
  )


### detecter les variants

vcf_long <- vcf_long %>%
  mutate(Variant_present =
           ifelse(Genotype == "./.:.", 0, 1))

### creer le passsage

vcf_long <- vcf_long %>%
  mutate(Passage = sub("-.*", "", Sample))

# Compter variants par passage
variants_par_passage <- vcf_long %>%
  filter(Variant_present == 1) %>%
  group_by(Passage, Sample) %>%
  summarise(Nb_variants = n(),
            .groups = "drop")


# Compter variants par passage

variants_par_passage_stat <- variants_par_passage %>%
  group_by(Passage) %>%
  summarise(
    mediane = median(Nb_variants, na.rm = TRUE),
    Q1 = quantile(Nb_variants, 0.25, na.rm = TRUE),
    Q3 = quantile(Nb_variants, 0.75, na.rm = TRUE),
    minimum = min(Nb_variants, na.rm = TRUE),
    maximum = max(Nb_variants, na.rm = TRUE),
    .groups = "drop"
  )




head(variants_par_passage)
str(variants_par_passage)

class(variants_par_passage$Nb_variants)

############# serie temporelle

# Transformer en facteur si ce n'est pas déjà fait
variants_par_passage$Sample <- factor(variants_par_passage$Sample, 
                                      levels = unique(variants_par_passage$Sample))  # garde l'ordre original
n_levels <- length(levels(variants_par_passage$Sample))
# Vérifie n_levels
print(n_levels)

# Créer un index numérique pour l'axe X
variants_par_passage$Index <- 1:nrow(variants_par_passage)

vline_pos <- which(levels(variants_par_passage$Sample) == "P25-5")
vline_pos2 <- which(levels(variants_par_passage$Sample) == "P27-5")


# Supposons que ta colonne Sample est déjà factor avec l'ordre correct
# On choisit de ne garder qu'un label sur deux
ggplot(variants_par_passage, aes(x = Sample, y = Nb_variants, group=1)) +
  geom_line(size = 1.2) +
  #geom_point(size = 3) +
  geom_vline(xintercept = c(vline_pos,vline_pos2), linetype = "dashed", color = "green", size = 1) +
  theme_classic() +
  scale_x_discrete(breaks = levels(variants_par_passage$Sample)[seq(1, length(levels(variants_par_passage$Sample)), by=2)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Evolution des variants par passage: Choc Thermique Froid",
       x = "Passages",
       y = "Nombre de variants")


################ boxplot des passages

ggplot(variants_par_passage, aes(x = Passage, y = Nb_variants)) +
  geom_boxplot() +
  theme_classic() +
  labs(title="",
       x = "Passages",
       y = "Nombre de variants")

ggsave(filename = "graphefroid.png", width = 8, height = 6, dpi = 300)


kruskal.test(Nb_variants ~ Passage, data = variants_par_passage)
pairwise.wilcox.test(
  variants_par_passage$Nb_variants,
  variants_par_passage$Passage,
  p.adjust.method = "BH"
)







########## Analyse des variants specifique, partagé 

#### creation de variable variants

# Exemple si ton dataframe s'appelle variants
vcf_long$Variant <- paste(vcf_long$POS, vcf_long$REF, vcf_long$ALT, sep = "-")


variant_summary <- vcf_long %>%
  filter(Variant_present == 1) %>%                  # on ignore les absents
  group_by(Variant) %>%
  summarise(
    n_passages = n_distinct(Passage),                  # combien de passages
    passages = paste(unique(Passage), collapse = ", ") # liste des passages
  ) %>%
  mutate(type = case_when(
    n_passages == 1 ~ "Spécifique",
    n_passages > 1  ~ "Partagé"
  ))


# Joindre le type à ton dataframe original
vcf_long <- vcf_long %>%
  left_join(variant_summary, by = "Variant")
table(vcf_long$type)

## Barplot des variants spécifiques vs partagés
ggplot(vcf_long %>% distinct(Variant, type), aes(x = type)) +
  geom_bar(fill = c("red", "blue")) +
  theme_classic() +
  labs(title = "Variants spécifiques vs partagés: froid", y = "Nombre de variants", x = "")


####### Visualiser
# Créer une matrice binaire Variant x Passage

variants_matrix <- vcf_long %>%
  select(Variant, Sample, Variant_present) %>%  # ne garder que ce qui sert
  pivot_wider(
    names_from = Sample, 
    values_from = Variant_present, 
    values_fill = 0
  )
str(variants_matrix)

variants_substress <- variants_matrix[, 
                                      colnames(variants_matrix) == "Variant" | 
                                        grepl("^(P25|P27)", colnames(variants_matrix))
]

head(variants_substress)

library(dplyr)

variants_global <- variants_substress %>%
  mutate(
    P25_global = ifelse(rowSums(select(., starts_with("P25"))) > 0, 1, 0),
    P27_global = ifelse(rowSums(select(., starts_with("P27"))) > 0, 1, 0)
  )

head(variants_global)

### extraire les catégorie biologiquement important

variants_partages <- variants_global %>%
  filter(P25_global == 1 & P27_global == 1)

head(variants_partages)

##### Spécifiques P25
variants_spec_P25 <- variants_global %>%
  filter(P25_global == 1 & P27_global == 0)

head(variants_spec_P25)
##### Spécifiques P27
variants_spec_P27 <- variants_global %>%
  filter(P25_global == 0 & P27_global == 1)
head(variants_spec_P27)

#### le nombre
nrow(variants_partages)
nrow(variants_spec_P25)
nrow(variants_spec_P27)

library(VennDiagram)
library(grid)

venn_list <- list(
  P25 = variants_spec_P25$Variant %>% union(variants_partages$Variant),
  P27 = variants_spec_P27$Variant %>% union(variants_partages$Variant)
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1.5
)

grid::grid.draw(venn.plot)

ggsave(filename = "partage_froid.png", width = 8, height = 6, dpi = 300)


library(UpSetR)

upset(fromList(venn_list),
      order.by = "freq",
      main.bar.color = "#E64B35",
      sets.bar.color = "#4DBBD5",
      matrix.color = "black",
      text.scale = 1.3,
      mainbar.y.label = "Nombre de variants",
      sets.x.label = "Taille des ensembles")

grid::grid.text(
  "Variants partagés et spécifiques entre les passages P25 et P27: Froid",
  x = 0.5, y = 0.98,
  gp = grid::gpar(fontsize = 16, fontface = "bold")
)


################ dans l'ensemble des base
variants_global_all <- variants_matrix %>%
  mutate(
    P15 = ifelse(rowSums(select(., starts_with("P15"))) > 0, 1, 0),
    P25 = ifelse(rowSums(select(., starts_with("P25"))) > 0, 1, 0),
    P27 = ifelse(rowSums(select(., starts_with("P27"))) > 0, 1, 0),
    P30 = ifelse(rowSums(select(., starts_with("P30"))) > 0, 1, 0),
    P50 = ifelse(rowSums(select(., starts_with("P50"))) > 0, 1, 0),
    P65 = ifelse(rowSums(select(., starts_with("P65"))) > 0, 1, 0),
    P90 = ifelse(rowSums(select(., starts_with("P90"))) > 0, 1, 0)
  )



names(variants_global_all)

variants_global_all <- variants_global_all[,c(1,37:43)]

venn_list_all <- list()

for (sample in colnames(variants_global_all)[-1]) {
  venn_list_all[[sample]] <- variants_global_all$Variant[variants_global_all[[sample]] == 1]
}

head(venn_list_all)

upset(fromList(venn_list_all), order.by = "freq",main.bar.color = "red")

library(UpSetR)
library(grid)

upset(fromList(venn_list_all),
      order.by = "freq",
      main.bar.color = "#E64B35",
      sets.bar.color = "#4DBBD5",
      matrix.color = "black",
      text.scale = 1.3,
      mainbar.y.label = "Nombre de variants",
      sets.x.label = "Taille des ensembles")

grid::grid.text(
  "Variants partagés et spécifiques entre les passages: Froid",
  x = 0.5, y = 0.98,
  gp = grid::gpar(fontsize = 16, fontface = "bold")
)


######### clustering

# Supposons que variants_matrix = lignes = variants, colonnes = Sample
# 0/1 ou fréquence

# Mettre la colonne Variant en rownames et garder seulement les colonnes numériques
# Suppose variants_matrix est un tibble wide
variants_matrix <- as.data.frame(variants_global_all)  # convertir tibble -> data.frame
rownames(variants_matrix) <- variants_matrix$Variant
variants_matrixclust <- as.matrix(variants_matrix[, -1])

# Assurer que c’est bien numérique
variants_matrixclust <- apply(variants_matrixclust, 2, as.numeric)

####### ACP

# Enlever la colonne Variant
mat <- variants_matrix[, -1]

# Mettre les variants comme noms de lignes
rownames(mat) <- variants_matrix$Variant

# Transposer → samples deviennent lignes
mat_t <- t(mat)


head (mat)

res_pca <- prcomp(mat_t, scale. = TRUE)



# Calcul de distance (ici distance euclidienne)
dist_variants <- dist(variants_matrixclust)

# Clustering hiérarchique
hc <- hclust(dist_variants, method = "ward.D2")

# Visualiser le dendrogramme
plot(hc, labels = rownames(variants_matrix), main = "Clustering des variants")
library(pheatmap)

variants_matrix$Variant <- NULL
variants_matrix <- as.matrix(variants_matrix)
variants_matrix <- apply(variants_matrix, 2, as.numeric)



pheatmap(variants_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white","red"))(50),
         show_rownames = FALSE,  # ou TRUE si peu de variants
         main = "Heatmap + Clustering des variants")




######## phylogenie

head(variants_matrix)
str(variants_matrix)


#### Calculer une distance entre échantillons
# variants_matrix = lignes = variants, colonnes = passages
# Conversion en numérique si besoin
# Convertir toutes les colonnes en numérique (présence/absence ou fréquence)
# Si variants_matrix a les rownames déjà

# Calcul de la distance (ex : distance binaire pour présence/absence)
dist_samples <- dist(t(variants_matrix), method = "binary")  # transposer pour distance entre colonnes


hc_samples <- hclust(dist_samples, method = "average")  # UPGMA
plot(hc_samples, main = "Phylogénie des passages selon les variants")

# Clustering hiérarchique
hc <- hclust(dist_variants, method = "ward.D2")

# Visualiser le dendrogramme
plot(hc, labels = rownames(hc_samples), main = "Clustering des variants")


# Clustering
hc_samples <- hclust(dist_samples, method = "ward.D2")

# Plot propre
plot(hc_samples, cex = 1.2, main = "Clustering des passages")


library(ape)

# Transformer en objet phylo
phylo_tree <- as.phylo(hc_samples)

# Dessiner l’arbre
plot(phylo_tree, main = "Phylogénie des passages: froid")

########## evolution uniquement entre les echantilllon de P25 et P27
# Vérifier les noms des colonnes
colnames(variants_matrix)

# Garder uniquement P25 et P27
variants_sub <- variants_matrix[, grepl("^P25|^P27", colnames(variants_matrix))]

#### verification
variants_sub <- as.data.frame(variants_sub)
variants_sub[] <- lapply(variants_sub, function(x) as.numeric(as.character(x)))
variants_sub[is.na(variants_sub)] <- 0

variants_sub <- as.matrix(variants_sub)

##### calculer la distace

dist_samples <- dist(t(variants_sub), method = "binary")



hc <- hclust(dist_samples, method = "average")  # UPGMA
tree <- as.phylo(hc)

plot(tree,
     main = "Phylogénie P25 vs P27: froid",
          cex = 0.8)

### ACP 

library(factoextra)

variants_sub


# Mettre les variants comme rownames si besoin
rownames(variants_sub) <- variants_sub$Variant

# Transposer → échantillons deviennent lignes
mat_t <- t(mat)


group <- ifelse(grepl("^P25", rownames(mat_t)), "P25", "P27")

fviz_pca_ind(res_pca,
             geom.ind = "point",
             col.ind = group,
             palette = c("blue","red"),
             addEllipses = TRUE,
             legend.title = "Passage",
             title = "ACP des variants entre P25 et P27")





####### calcule de distance SNP
## Si tu veux une distance génétique plus spécifique :
library(dartR)
# snps_genlight : objet genlight avec les SNP
class(snps_genlight)

dist_snp <- utils.dist.ind.snp(snps_genlight, method = "Euclidean")  # ou autre


library(vegan)

# data.frame de méta‑données
meta <- data.frame(groupe = variants_sub)

# PERMANOVA
res <- adonis2(dist_snp ~ variants_sub, data = meta, permutations = 999)
print(res)

#### Pour « voir » la séparation P25/P27, tu peux faire une PCoA sur la même matrice de distance :
pcoa_res <- cmdscale(dist_snp, k = 2, eig = TRUE)
scores <- as.data.frame(pcoa_res$points)
scores$groupe <- groupe

library(ggplot2)
ggplot(scores, aes(Dim1, Dim2, color = groupe)) +
  geom_point(size = 3) +
  theme_minimal()




############################################################ Analyse chaud ###########
library(dplyr)
library(tidyr)
library(ggplot2)
library(UpSetR)
library(magrittr)
library(reshape2)
# Lire le TSV
vcf_chaud <- read.delim("Base_chaud.tsv", stringsAsFactors = FALSE, check.names = FALSE)

# Vérifier les colonnes
head(vcf_chaud)

#### DIFFERENCIÉ les snp ,insert, deletion

vcf_chaud$typevariant <-ifelse(nchar(vcf_chaud$REF) < nchar(vcf_chaud$ALT), "Insertion" , 
                         ifelse (nchar(vcf_chaud$REF) > nchar (vcf_chaud$ALT), "deletion",
                                 "SNP"))

table_variants <- table(vcf_chaud$typevariant)
table_variants

print(vcf_chaud[vcf_chaud$typevariant == "Insertion", ])

##### metre mes sample en colone
vcf_long <- vcf_chaud %>%
  pivot_longer(
    cols = -c(CHROM:FORMAT, typevariant),
    names_to = "Sample",
    values_to = "Genotype"
  )


### detecter les variants

vcf_long <- vcf_long %>%
  mutate(Variant_present =
           ifelse(Genotype == "./.:.", 0, 1))

### creer le passsage

vcf_long <- vcf_long %>%
  mutate(Passage = sub("-.*", "", Sample))

# Compter variants par passage
variants_par_passage <- vcf_long %>%
  filter(Variant_present == 1) %>%
  group_by(Passage, Sample) %>%
  summarise(Nb_variants = n(),
            .groups = "drop")



variants_par_passage$Sample


head(variants_par_passage)
str(variants_par_passage)

class(variants_par_passage$Nb_variants)

############# serie temporelle



variants_par_passage <- variants_par_passage %>%
  # extraire le numéro après le tiret
  mutate(num = as.numeric(sub(".*-", "", Sample))) %>%
  # trier d'abord par Passage, puis par numéro
  arrange(Passage, num) %>%
  select(-num)  # enlever la colonne temporaire si pas besoin

variants_par_passage$Sample
# Transformer en facteur si ce n'est pas déjà fait
variants_par_passage$Sample <- factor(variants_par_passage$Sample, 
                                      levels = unique(variants_par_passage$Sample))  # garde l'ordre original
n_levels <- length(levels(variants_par_passage$Sample))
# Vérifie n_levels
print(n_levels)

# Créer un index numérique pour l'axe X
variants_par_passage$Index <- 1:nrow(variants_par_passage)

# Trouver la position numérique de P25-10
vline_pos <- which(levels(variants_par_passage$Sample) == "P25-10")
vline_pos2 <- which(levels(variants_par_passage$Sample) == "P27-10")

# Supposons que ta colonne Sample est déjà factor avec l'ordre correct
# On choisit de ne garder qu'un label sur deux
ggplot(variants_par_passage, aes(x = Sample, y = Nb_variants, group=1)) +
  geom_line(size = 1.2) +
  #geom_point(size = 3) +
  geom_vline(xintercept = c(vline_pos,vline_pos2), linetype = "dashed", color = "red", size = 1) +
  theme_classic() +
  scale_x_discrete(breaks = levels(variants_par_passage$Sample)[seq(1, length(levels(variants_par_passage$Sample)), by=2)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Evolution des variants par passage: Choc Thermique chaud",
       x = "Passages",
       y = "Nombre de variants")




################ boxplot des passages

ggplot(variants_par_passage, aes(x = Passage, y = Nb_variants)) +
  geom_boxplot() +
  theme_classic() +
  labs(title="",
       x = "Passages",
       y = "Nombre de variants")

ggsave(filename = "boxplot.png", width = 8, height = 6, dpi = 300)


########## Analyse des variants specifique, partagé 

#### creation de variable variants

# Exemple si ton dataframe s'appelle variants
vcf_long$Variant <- paste(vcf_long$POS, vcf_long$REF, vcf_long$ALT, sep = "-")


variant_summary <- vcf_long %>%
  filter(Variant_present == 1) %>%                  # on ignore les absents
  group_by(Variant) %>%
  summarise(
    n_passages = n_distinct(Passage),                  # combien de passages
    passages = paste(unique(Passage), collapse = ", ") # liste des passages
  ) %>%
  mutate(type = case_when(
    n_passages == 1 ~ "Spécifique",
    n_passages > 1  ~ "Partagé"
  ))


# Joindre le type à ton dataframe original
vcf_long <- vcf_long %>%
  left_join(variant_summary, by = "Variant")

## Barplot des variants spécifiques vs partagés
ggplot(vcf_long %>% distinct(Variant, type), aes(x = type)) +
  geom_bar(fill = c("red", "blue")) +
  theme_classic() +
  labs(title = "Variants spécifiques vs partagés:chaud", y = "Nombre de variants", x = "")

####### Visualiser
# Créer une matrice binaire Variant x Passage

variants_matrix <- vcf_long %>%
  select(Variant, Sample, Variant_present) %>%  # ne garder que ce qui sert
  pivot_wider(
    names_from = Sample, 
    values_from = Variant_present, 
    values_fill = 0
  )
str(variants_matrix)

variants_substress <- variants_matrix[, 
                                      colnames(variants_matrix) == "Variant" | 
                                        grepl("^(P25|P27)", colnames(variants_matrix))
]

head(variants_substress)

library(dplyr)

variants_global <- variants_substress %>%
  mutate(
    P25_global = ifelse(rowSums(select(., starts_with("P25"))) > 0, 1, 0),
    P27_global = ifelse(rowSums(select(., starts_with("P27"))) > 0, 1, 0)
  )

head(variants_global)

### extraire les catégorie biologiquement important

variants_partages <- variants_global %>%
  filter(P25_global == 1 & P27_global == 1)

head(variants_partages)

##### Spécifiques P25
variants_spec_P25 <- variants_global %>%
  filter(P25_global == 1 & P27_global == 0)

head(variants_spec_P25)
##### Spécifiques P27
variants_spec_P27 <- variants_global %>%
  filter(P25_global == 0 & P27_global == 1)
head(variants_spec_P27)

#### le nombre
nrow(variants_partages)
nrow(variants_spec_P25)
nrow(variants_spec_P27)

library(VennDiagram)

venn_list <- list(
  P25 = variants_spec_P25$Variant %>% union(variants_partages$Variant),
  P27 = variants_spec_P27$Variant %>% union(variants_partages$Variant)
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 1.5
)

grid::grid.draw(venn.plot)

ggsave(filename = "graph2chaud.png", width = 8, height = 6, dpi = 300)


library(UpSetR)


upset(fromList(venn_list),
      order.by = "freq",
      main.bar.color = "#E64B35",
      sets.bar.color = "#4DBBD5",
      matrix.color = "black",
      text.scale = 1.3,
      mainbar.y.label = "Nombre de variants",
      sets.x.label = "Taille des ensembles")

grid::grid.text(
  "Variants partagés et spécifiques entre les passages P25 et P27: chaud",
  x = 0.5, y = 0.98,
  gp = grid::gpar(fontsize = 16, fontface = "bold")
)


################ dans l'ensemble des base
variants_global_all <- variants_matrix %>%
  mutate(
    P15_global = ifelse(rowSums(select(., starts_with("P15"))) > 0, 1, 0),
    P25_global = ifelse(rowSums(select(., starts_with("P25"))) > 0, 1, 0),
    P27_global = ifelse(rowSums(select(., starts_with("P27"))) > 0, 1, 0),
    P30_global = ifelse(rowSums(select(., starts_with("P30"))) > 0, 1, 0),
    P50_global = ifelse(rowSums(select(., starts_with("P50"))) > 0, 1, 0),
    P65_global = ifelse(rowSums(select(., starts_with("P65"))) > 0, 1, 0),
    P90_global = ifelse(rowSums(select(., starts_with("P90"))) > 0, 1, 0)
  )



names(variants_global_all)

variants_global_all <- variants_global_all[,c(1,37:43)]

venn_list_all <- list()

for (sample in colnames(variants_global_all)[-1]) {
  venn_list_all[[sample]] <- variants_global_all$Variant[variants_global_all[[sample]] == 1]
}

head(venn_list_all)



library(UpSetR)
library(grid)

upset(fromList(venn_list_all),
      order.by = "freq",
      main.bar.color = "#E64B35",
      sets.bar.color = "#4DBBD5",
      matrix.color = "black",
      text.scale = 1.3,
      mainbar.y.label = "Nombre de variants",
      sets.x.label = "Taille des ensembles")

grid::grid.text(
  "Variants partagés et spécifiques entre les passages: Froid",
  x = 0.5, y = 0.98,
  gp = grid::gpar(fontsize = 16, fontface = "bold")
)



######### clustering

# Supposons que variants_matrix = lignes = variants, colonnes = Sample
# 0/1 ou fréquence

# Mettre la colonne Variant en rownames et garder seulement les colonnes numériques
# Suppose variants_matrix est un tibble wide
variants_matrix <- as.data.frame(variants_matrix)  # convertir tibble -> data.frame
rownames(variants_matrix) <- variants_matrix$Variant
variants_matrixclust <- as.matrix(variants_matrix[, -1])

# Assurer que c’est bien numérique
variants_matrixclust <- apply(variants_matrixclust, 2, as.numeric)


# Calcul de distance (ici distance euclidienne)
dist_variants <- dist(variants_matrixclust)

# Clustering hiérarchique
hc <- hclust(dist_variants, method = "ward.D2")

# Visualiser le dendrogramme
plot(hc, labels = rownames(variants_matrix), main = "Clustering des variants")
library(pheatmap)

variants_matrix$Variant <- NULL
variants_matrix <- as.matrix(variants_matrix)
variants_matrix <- apply(variants_matrix, 2, as.numeric)

library(pheatmap)

pheatmap(variants_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white","red"))(50),
         show_rownames = FALSE,  # ou TRUE si peu de variants
         main = "Heatmap + Clustering des variants")


######## phylogenie

head(variants_matrix)
str(variants_matrix)


#### Calculer une distance entre échantillons
# variants_matrix = lignes = variants, colonnes = passages
# Conversion en numérique si besoin
# Convertir toutes les colonnes en numérique (présence/absence ou fréquence)
# Si variants_matrix a les rownames déjà

# Calcul de la distance (ex : distance binaire pour présence/absence)
dist_samples <- dist(t(variants_matrix), method = "binary")  # transposer pour distance entre colonnes


hc_samples <- hclust(dist_samples, method = "average")  # UPGMA
plot(hc_samples, main = "Phylogénie des passages selon les variants")


library(ape)

# Transformer en objet phylo
phylo_tree <- as.phylo(hc_samples)

# Dessiner l’arbre
plot(phylo_tree, main = "Phylogénie des passages: chaud")

########## evolution uniquement entre les echantilllon de P25 et P27
# Vérifier les noms des colonnes
colnames(variants_matrix)

# Garder uniquement P25 et P27
variants_sub <- variants_matrix[, grepl("^P25|^P27", colnames(variants_matrix))]

#### verification
variants_sub <- as.data.frame(variants_sub)
variants_sub[] <- lapply(variants_sub, function(x) as.numeric(as.character(x)))
variants_sub[is.na(variants_sub)] <- 0

variants_sub <- as.matrix(variants_sub)

##### calculer la distace

dist_samples <- dist(t(variants_sub), method = "binary")



hc <- hclust(dist_samples, method = "average")  # UPGMA
tree <- as.phylo(hc)

plot(tree,
     main = "Phylogénie P25 vs P27: chaud",
     cex = 0.8)

  
  









  




