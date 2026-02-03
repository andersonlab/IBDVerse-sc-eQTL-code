library(tidyverse)
library(ggrepel)
library(ggpattern)
library(scatterpie)
library(patchwork)

##############
# Variables
###############

server <- 'farm'
# server <- 'openstack'
tissue <- 'multi_tissue'


if (server == 'farm'){
  repo.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/'
  interaction.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_12-multi_tissue_interaction_results/TensorQTL_eQTLS/'
  out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/'
  vcf.path <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/proc_data/genotypes/Aug2024_ashg/CCF_OTAR-plates_12345_imputed_allchr-no_rsid_expr_filt.vcf.gz'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  sumstats.all.basedir <- paste0('/home/rstudio/eqtl/data/IBDverse/',tissue,'/pseudobulk-base/')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/coloc/'
  
}

data.dir <- paste0(repo.dir,'/data/')
normalised <- FALSE

# Plotting
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE

##############
# CODE
##############

interaction.sumstat.df <- read_ieqtls(interaction.basedir)

sig.interaction.sumstat.df <- interaction.sumstat.df %>% 
  filter(pval_adj_bh < 0.05,
         tissue != 'Cross-site')

# Load in the interaction variable dataframe so I can extract disease status
# for each individual
interaction.variable.df <- read_tsv('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/interaction_files/2024_12_27-multi_tissue/eQTL_interactions.tsv')

# Read the IBD GWAS snps for aligning purposes
gwas.snps <- read_ibd_gwas('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/gwas_sumstats_final/')

# Make gene level box plot
interaction.effect <- 'ses_cd_binned'
interaction.labels <- c('0'='Uninflamed', '1'='Mild', '2'='Moderate/Severe')
gene.of.interest <- 'RNF14' # Enterocyte

eqtl.of.interest <- filter(sig.interaction.sumstat.df,
                           symbol == gene.of.interest,
                           interaction == interaction.effect)
variant.of.interest <- eqtl.of.interest$variant_id

cur.gwas.snp <- gwas.snps %>% 
  filter(RSid == variant.of.interest) %>% 
  slice_min(pval)


# Load expression data
if (!normalised){
  raw.expression.paths <- list.files(path = paste0(sumstats.all.basedir, "../aggregated_counts/"), 
                                 pattern = "phenotype_file.tsv$",
                                 recursive = TRUE,  full.names = TRUE)
  expression.df <- read.table(str_subset(raw.expression.paths, eqtl.of.interest$annotation),header = TRUE) %>% 
    rename(ensgene = ENS)
} else {
  int.expression.paths <- list.files(path = paste0(sumstats.all.basedir, "../norm_data/"), 
                                 pattern = "normalised_phenotype.tsv$",
                                 recursive = TRUE,  full.names = TRUE)
  expression.df <- read.table(str_subset(int.expression.paths, eqtl.of.interest$annotation), header = TRUE, row.names = NULL) %>% 
    rename(ensgene = row.names)
}



# Load vcf
if (!exists('vcf.df')){
  print('Loading VCF')
  vcf.df <- read.table(vcf.path,skip=52, header=TRUE, comment.char = "")
}


variant.df <- filter(vcf.df, ID == variant.of.interest) %>% 
  pivot_longer(!names(vcf.df)[1:9] , names_to = 'genotype_id') %>% 
  mutate(genotype_id = gsub('^X', '', genotype_id)) %>% 
  separate(value, into = str_split('GT:DS:HDS:GP', ':')[[1]],sep = ':') %>% 
  mutate(Genotype = str_count(GT, "1")) %>% 
  select(ID,genotype_id,Genotype,GT) %>% 
  left_join(interaction.variable.df, by = c('genotype_id'='Genotyping_ID')) %>% 
  drop_na() 


# Slice for a particular gene
gene.df <- filter(expression.df,ensgene == eqtl.of.interest$phenotype_id) %>% 
  pivot_longer(!ensgene, names_to = 'annotation_individual') %>% 
  separate(annotation_individual,
           into = c('annotation_type','annotation', 'individual') ,sep = '[.]') %>% 
  separate(individual, c('method','genotype_id', 'sanger_sample_id'), sep='_',extra='merge') %>% 
  left_join(variant.df, by='genotype_id') %>% 
  drop_na() # Why do I need this, check data

if (cur.gwas.snp$beta < 0){
  gene.df <- gene.df %>%
    mutate(Genotype = abs(Genotype - 2))
}


# Genotype by Expression plot
ggplot(gene.df, aes(y=value, x=Genotype, colour=as.factor(ses_cd_binned),
                    fill=as.factor(Genotype))) +
  geom_boxplot(outlier.colour = NA) +
  geom_abline(slope = 0, intercept = 0, size = 0.2) +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             size=2, shape=19,stroke=1, alpha=0.6) +
  scale_fill_grey(start = 1, end = 1,guide = 'none') +
  scale_colour_viridis_d(option='plasma', end = 0.7,
                         labels = interaction.labels) +
  scale_x_continuous(breaks = seq(0, 2, 1), 
                     labels = c('0'='0', '1'='1', '2'='2')) +
  theme_classic() + 
  labs(title = paste0(eqtl.of.interest$label_new,'\n', eqtl.of.interest$variant_id),
       x=paste('Number of',cur.gwas.snp$Disease, 'risk increasing alleles'),
       y = paste0("Pseudobulked expression of\n",eqtl.of.interest$symbol),
       color = 'Inflammation status')

if (SAVE.PLOTS){
  ggsave(paste0(out.dir,'interaction-',eqtl.of.interest$interaction_new,'-',
                eqtl.of.interest$label_new,'-',
                eqtl.of.interest$symbol,
                '_genotype_boxplot.pdf'),
         device = cairo_pdf,width=7, height=4)
}

# Interaction by expression plot
ggplot(gene.df, aes(y=value, x=as.factor(ses_cd_binned), fill=as.factor(Genotype),
                    colour=as.factor(ses_cd_binned))) +
  geom_abline(slope = 0, intercept = 0, size = 0.2,linetype='dashed') +
  geom_point(position=position_jitterdodge(dodge.width=0.75),
             size=3, shape=19,stroke=1, alpha=0.6) +
  geom_boxplot(outlier.colour = NA, alpha=0.6, lwd=1) +
  scale_fill_grey(start = 0.9, end = 0.1) +
  scale_colour_viridis_d(option='plasma', end = 0.7,
                         labels = interaction.labels) +
  theme_classic() +
  scale_x_discrete(labels=interaction.labels) +
  labs(title = paste0(eqtl.of.interest$label_new,'\n',
                      eqtl.of.interest$variant_id,'\n',
                      'AF = ', round(eqtl.of.interest$af, 3), '   ',
                      'adj.P = ', round(eqtl.of.interest$pval_adj_bh, 4)),
       x = interaction.effect,
       y = paste0("Normalised expression of\n",eqtl.of.interest$symbol),
       fill = "Minor allele\ncount",
       colour = gsub(' ', '\n',interaction.effect)) +
  theme(title = element_text(size=13),
        axis.title = element_text(size=12),
        legend.position = 'bottom'
  ) +
  guides(color = guide_legend(order = 2),fill = guide_legend(order = 1))
