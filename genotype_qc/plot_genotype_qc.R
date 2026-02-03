library(tidyverse)
library(ggrepel)
library(patchwork)


# working.dir <- "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/genotypes/dna_genotypes/2023_Nov/"
working.dir <- "/home/rstudio//eqtl/data/genotypes_2024_march/"
plots.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/genetics/'
SAVE.PLOTS <- TRUE
setwd(working.dir)

id.mapping <- read_table('plate_1-5_mapping.txt', col_names = c('original', 'new'))

gut.metadata <- read_csv('GUT_scRNAseq_metadata - GUT_scRNAseq-cleaned.csv')

genotype_plate_batch_effects <- gut.metadata %>% 
  dplyr::select(Pre_QC_Genotyping_ID, Genotyping_material, Genotyping_status) %>% 
  distinct() %>% 
  filter(str_detect(Genotyping_status,'^Batch')) %>% 
  mutate(Pre_QC_Genotyping_ID = str_replace(Pre_QC_Genotyping_ID,'Plate ', 'Batch3_')) %>%
  add_row(Pre_QC_Genotyping_ID = '3981576863664',  # Adding two samples that were mixed up manually
          Genotyping_status = 'Batch 1, sent Dec 2021',
          Genotyping_material = 'Blood') %>% 
  add_row(Pre_QC_Genotyping_ID = '3981576864678',
          Genotyping_status = 'Batch 1, sent Dec 2021',
          Genotyping_material = 'Blood') %>% 
  mutate(Genotyping_material = str_replace(Genotyping_material,
                                           'Tissue[/]Cell suspension',
                                           'Cell Suspension'))   # Remove the word tissue 

# Principal component analysis

#Read in eigenvecs
pca_eigenvecs = read.table("pca_pruned.eigenvec", header = F) %>% 
  mutate(genotype_id = recode(V1, !!!deframe(id.mapping))) %>% 
  left_join(genotype_plate_batch_effects, by=c('genotype_id'='Pre_QC_Genotyping_ID')) %>% 
  arrange(desc(Genotyping_material))

# Eigenvalues
pca_eigenvals = read.table('pca_pruned.eigenval', header=F) %>% 
  mutate(max_var_explained = V1/sum(V1))

# Make the plots (dropping harriet samples)
p1 <- ggplot(drop_na(pca_eigenvecs), aes(x=V3, y=V4, label = genotype_id,
                          colour=Genotyping_status, shape = Genotyping_material)) + 
  geom_point(size=4, alpha=0.4,) +
  # geom_text_repel() + 
  theme_classic() + 
  labs(x = 'Within-study PC1', y = 'Within-study PC2') +
  theme(legend.position='none')

p2 <- ggplot(drop_na(pca_eigenvecs), aes(x=V5, y=V6, label = genotype_id,
                          colour=Genotyping_status, shape = Genotyping_material)) + 
  geom_point(size=4, alpha=0.4) +
  # geom_text_repel() + 
  theme_classic() + 
  labs(x = 'Within-study PC3', y = 'Within-study PC4')

# Patch together the plots
p1 + p2

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(plots.dir,"within_sample_genetic_pca.pdf"),
         width = 9, height = 4,device = cairo_pdf)
}

#### Project into 1000G PC space and predict ancestry #####

projected_data <- read.table(file="1000g_projectedpc.txt", header=T) %>%
  mutate(genotype_id = recode(FID, !!!deframe(id.mapping))) %>% 
  left_join(read.table('1000g_projected_popref.txt', header=TRUE)) %>% 
  left_join(read.table(file='1000g_projected_InferredAncestry.txt',
                              header=T)) %>% 
  mutate(Study = AFF == 2) %>% 
  dplyr::select(genotype_id, FID, PC1, PC2, PC3, PC4, Study, Ancestry, Population, Pr_1st)

ancestry.color.palette <- c("EUR" = "#1f77b4",
                            "EAS" = "#ff7f0e",
                            "AMR" = "#2ca02c",
                            "SAS" = "#d62728",
                            "AFR" = "#F0E442",
                            "AMR;EUR" = "#258b70",
                            "EUR;SAS" = "#7a4f6e")


# Plot samples overlaid on 1000G ancestries
pca.1000g.left <- ggplot(filter(projected_data, Study == FALSE),
                         aes(x=PC1, y=PC2, colour=Population)) +
  geom_point(alpha=0.5) + 
  theme_classic() +
  scale_color_manual(values=ancestry.color.palette,
                     guide = guide_legend(override.aes = list(size = 3,
                                                              alpha = 1) )) +
  theme(legend.position = "none")+
  labs(x = '1000G PC1', y = '1000G PC2')

pca.1000g.right <- ggplot(filter(projected_data, Study == FALSE),
                          aes(x=PC3, y=PC4, colour=Population)) +
  geom_point(alpha=0.5) + 
  theme_classic() +
  scale_color_manual(values=ancestry.color.palette,
                     guide = guide_legend(override.aes = list(size = 3,
                                                              alpha = 1) )) +
  theme(legend.position = "right")+
  labs(x = '1000G PC3', y = '1000G PC4')

pca.1000g.left + pca.1000g.right

# Plot predicted ancestries
projected.pca.left <- ggplot(filter(projected_data, Study == TRUE),
                             aes(x=PC1, y=PC2, colour=Ancestry)) +
  geom_point(data=projected_data,alpha=0) +
  geom_point(alpha=0.75) + 
  geom_text_repel(data = filter(projected_data, Ancestry != 'EUR'),
                  aes(label = Ancestry),show.legend=F) +
  theme_classic() +
  guides(alpha='none') +
  scale_color_manual(values=ancestry.color.palette,
                     guide = guide_legend(override.aes = list(size = 3,
                                                              alpha = 1)),
                     name = "Predicted\npopulation") +
  theme(legend.position = "none")+
  labs(x = '1000G PC1', y = '1000G PC2')

projected.pca.right <- ggplot(filter(projected_data, Study == TRUE),
                              aes(x=PC3, y=PC4)) +
  geom_point(data=projected_data, alpha=0) +
  geom_point(aes(colour=Ancestry), alpha=0.75) + 
  geom_text_repel(data = filter(projected_data, Study == TRUE, Ancestry != 'EUR'),
                  aes(label = Ancestry,colour=Ancestry),show.legend=F) +
  theme_classic() +
  guides(alpha='none') +
  scale_color_manual(values=ancestry.color.palette,
                     guide = guide_legend(override.aes = list(size = 3,
                                                              alpha = 1)),
                     name = "Predicted\npopulation") +
  theme(legend.position = "right") +
  labs(x = '1000G PC3', y = '1000G PC4')

(pca.1000g.left + pca.1000g.right) / (projected.pca.left + projected.pca.right) + 
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(plots.dir,"1000G_PCA_projection.pdf"),
         width = 8, height = 8,device = cairo_pdf)
}


## Plot missingness

sample_pre_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.imiss",head=T)
var_pre_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.lmiss",head=T)

sample_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.imiss",head=T)
var_miss<-read.table("merged_plates_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss",head=T)


ggplot(sample_miss, aes(x=F_MISS)) +
  geom_histogram() + 
  theme_classic() +
  labs(y='Number of samples', x='Missingness rate')

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(plots.dir,"sample_missingness.pdf"),
         width = 4, height = 4,device = cairo_pdf)
}

## Plot KING relatedness

kinship<-read.table("kinship.kin0",head=T) %>% 
  group_by(FID1) %>% 
  filter(Kinship == max(Kinship)) %>% 
  mutate(FID1 = recode(FID1, !!!deframe(id.mapping)),
         FID2 = recode(FID2, !!!deframe(id.mapping))) %>% 
  dplyr::select(-ID1, -ID2)

ggplot(kinship, aes(x=Kinship)) +
  geom_histogram() + 
  theme_classic() +
  labs(y='Number of samples', x='Kinship coefficient')

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(plots.dir,"kinship.pdf"),
         width = 4, height = 4,device = cairo_pdf)
}
