library(tidyverse)
library(patchwork)
library(ggh4x)

##############
# Variables
###############

server <- 'farm'
# server <- 'openstack'

if (server == 'farm'){
  repo.dir <- '~/eqtl/code/IBDVerse-sc-eQTL-code/'
  coloc.base.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed'
  coloc.interaction.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed'
  out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  coloc.dir <- paste0('/home/ubuntu/eqtl/data/IBDverse/multi_tissue/coloc-',interaction.or.base)
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/coloc/'
}


data.dir <- paste0(repo.dir,'data/')
otar.ibd.coloc.path <- paste0(data.dir,'all_otg_ibd_colocs.txt')
clumps.path <- paste0(data.dir,'clumped_all.txt.gz')
# SAVE
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE


##############
# Functions
##############

source(paste0(repo.dir,'/qtl_plot/helper_functions.R'))

##############
# Main
###############

# all.coloc.df <- get_colocs(coloc.dir)

known.base.coloc.df <- get_colocs(coloc.base.dir, known_ibd_only = TRUE)
known.interaction.coloc.df <- get_colocs(coloc.interaction.dir, known_ibd_only = TRUE) %>% 
  mutate(interaction = str_split_fixed(condition_name, "__", 3)[,3])

known.coloc.df <- bind_rows(mutate(known.base.coloc.df, interaction = 'base'),
                                 known.interaction.coloc.df)



top.coloc.hits <- known.coloc.df %>%
# top.coloc.hits <- all.coloc.df %>% mutate(loci_ID = gwas_lead) %>%
  filter(PP.H4.abf > PP4.loose.threshold) %>% 
  arrange(-PP.H4.abf) %>%
  group_by(gwas_trait) %>%
  distinct(loci_ID,.keep_all = TRUE) %>%
  ungroup() %>%
  arrange(-PP.H4.abf) %>% 
  mutate(category_new = factor(category_new, levels = sort(unique(category_new))))

otar.ibd.colocs <- read_otar(repo.dir)

plot.otar.top.hits <- otar.ibd.colocs %>%
  arrange(desc(H4)) %>% 
  distinct(loci_ID,.keep_all = TRUE) %>%
  distinct(Gene_symbol,.keep_all = TRUE)

otar.hits <- unique(filter(otar.ibd.colocs, H4 > PP4.loose.threshold)$Gene_symbol)

monogenic.ibd.genes <- read_monogenic(repo.dir) 

plot.top.coloc.hits <- top.coloc.hits %>% 
  mutate(in_OTAR_genetics = factor(phenotype_id %in% otar.hits,
                                   levels = c(TRUE,FALSE)),
         is_monogenic = factor(phenotype_id %in% monogenic.ibd.genes$ensgene,
                               levels = c(TRUE,FALSE))) %>% 
  filter(gwas_trait %in% ibd.traits) %>% 
  group_by(loci_ID) %>% 
  arrange(annotation_type) %>% 
  mutate(all_levels = paste(unique(annotation_type), collapse=', ')) %>% 
  ungroup() %>% 
  arrange(-PP.H4.abf) %>% 
  distinct(loci_ID,.keep_all = TRUE) %>% 
  distinct(phenotype_id,.keep_all = TRUE)

# Make a pie chart of the loci with OTAR colocs and the loci our colocs add
loci.321.pi.df <- tibble(loci_ID = plot.otar.top.hits$loci_ID, source='OTG',
                         label = 'Public QTL\ncolocalised loci') %>% 
  bind_rows(tibble(loci_ID = plot.top.coloc.hits$loci_ID, source='IBDverse', # Comment out these two
                   label = '\n\nIBDverse-only\ncolocalised loci')) %>% # Lines and re-run to "animate"
  bind_rows(tibble(loci_ID = ibd.regions$loci_ID, source='No colocs',
                   label = 'Remaining\nIBD loci')) %>%
  distinct(loci_ID,.keep_all = TRUE) %>% 
  add_count(source) %>% 
  mutate(source = factor(source, levels = c('No colocs','IBDverse','OTG')))#,
         # label = paste0(label,'\nn=',n))

loci.pi.p <- ggplot(loci.321.pi.df, aes(x = '', fill = source)) +
  geom_bar(colour='black') +
  geom_text(stat = 'count',aes(label = paste0(label,'\nn=',n),
                               y=after_stat(count),colour=source),
            position = position_stack(vjust = 0.5),
            size=2.2,fontface='bold',) +
  coord_polar("y", start = 0) + 
  # ylim(0, 321) +
  scale_y_continuous(breaks = cumsum(unique(loci.321.pi.df$n))) +
  scale_fill_manual(values=c('IBDverse' = otar.palette[["FALSE"]],
                             'OTG' =  otar.palette[["TRUE"]],
                             'No colocs' = 'white')) +
  scale_color_manual(values=c('IBDverse' = 'black',
                             'OTG' = 'black',
                             'No colocs' = 'black')) +
  guides(fill='none', colour='none') +
  theme_void() +
  theme(axis.text = element_text(size=8))

print(loci.pi.p)

# Scatterplot of relative contributions to proportion of eGenes and proportion of colocs
#NOTE: Needs df.count.sig.egenes from plot_n_egenes.R
egene.prop.df <- df.count.sig.egenes %>% 
  dplyr::count(annotation_type) %>%
  mutate(frequency = n / sum(n)) %>% 
  mutate(proportion_class = 'eGene')

# Coloc prop (if ANY coloc gene in locus can be found for tissue, cat. or cell type)
coloc.prop.df <- known.base.coloc.df %>% 
  filter(PP.H4.abf > PP4.loose.threshold,
         gwas_trait %in% ibd.traits) %>% 
  # filter(phenotype_id %in% plot.top.coloc.hits$phenotype_id) %>% 
  group_by(loci_ID) %>% 
  slice_min(Level,with_ties = FALSE) %>% 
  ungroup() %>% 
  dplyr::count(annotation_type) %>% 
  mutate(frequency = n / sum(n)) %>% 
  mutate(proportion_class = 'coloc')

prop.df <- bind_rows(egene.prop.df,coloc.prop.df) %>% 
  pivot_wider(names_from=proportion_class,values_from = c('frequency','n')) %>% 
  mutate(enrichment = round(frequency_coloc/frequency_eGene,2),
         annotation_type = factor(annotation_type,
                                  unname(level.mapping)))


p.coloc.enrichment <- ggplot(data=prop.df, aes(y=enrichment, x=frequency_eGene,fill=annotation_type,
                         label=paste0(enrichment,'x'))) +
  geom_abline(slope=0, intercept=1, linetype="dashed", colour='grey') +
  geom_segment(aes(y=1, xend=frequency_eGene, yend=enrichment)) +
  geom_point(size=5,colour='black',shape=21) +
  geom_text_repel(box.padding = 0.5, point.padding = 2, nudge_y = 0.08,
                  fontface = "bold", size=5,show.legend = FALSE,
                  min.segment.length = 3, colour='black',bg.colour='white') +
  scale_fill_manual(values=annot.class.palette, name="Granularity") +
  labs(x='Proportion of eGenes',
       y='Odds ratio') +
  theme_classic() +
  theme(legend.position = 'bottom')


# Now do a major population prop df
egene.annot.prop.df <- sig.sumstat.df %>% 
  filter(Level != 0) %>% 
  arrange(qval) %>%
  distinct(phenotype_id, .keep_all = TRUE) %>% 
  dplyr::count(category_new) %>%
  mutate(frequency = n / sum(n)) %>% 
  mutate(proportion_class = 'eGene')


# Coloc prop (if top coloc gene in locus can be found for tissue, cat. or cell type)
coloc.annot.prop.df <- known.base.coloc.df %>% 
  filter(PP.H4.abf > PP4.loose.threshold,
         gwas_trait %in% ibd.traits) %>% 
  filter(Level != 0) %>% 
  group_by(loci_ID) %>% 
  slice_min(PP.H4.abf,with_ties = FALSE) %>% 
  ungroup() %>% 
  dplyr::count(category_new) %>% 
  mutate(frequency = n / sum(n)) %>% 
  mutate(proportion_class = 'coloc')

annot.prop.df <- bind_rows(egene.annot.prop.df,coloc.annot.prop.df) %>% 
  pivot_wider(names_from=proportion_class,values_from = c('frequency','n')) %>% 
  mutate(enrichment = round(frequency_coloc/frequency_eGene,2))

ggplot(data=annot.prop.df, aes(y=enrichment, x=category_new,fill=category_new,
                         label=paste0(enrichment,'x'))) +
  geom_point(aes(size=n_coloc),colour='black',shape=21) +
  geom_text_repel(box.padding = 0.5, point.padding = 2, nudge_y = 0.08,
                  fontface = "bold", size=5,show.legend = FALSE,
                  min.segment.length = 3, colour='black',bg.colour='white') +
  geom_abline(slope=0, intercept=1, linetype="dashed") +
  scale_fill_manual(values=umap.category.palette, name="Granularity",
                    guide='none') + 
  scale_size(range=c(2, 10))+
  labs(y='Odds ratio') +
  guides(size = guide_legend("Number of\nIBD colocs"),) +
  theme_classic() +
  theme(legend.position = c(0.9,0.75),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=60,hjust=1))

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"coloc_eGene_annot_proportions.pdf"),
         width = 5, height = 3.7,device = cairo_pdf)
}


# Now read in the variance explained file so I can summarise which cell types
# are the likely effector cell types
if (!exists("var.explained.df")){
  var.explained.df <- read_tsv(paste0(out.dir,"colocs_table-var_explained-",
                                      gsub('[.]', 'pt', PP4.loose.threshold),".tsv"))
}

top.var.explained.df <- var.explained.df %>% 
  group_by(loci_ID) %>% 
  slice_max(rsquared, n=1) %>% 
  ungroup() %>% 
  annotate_sumstats() %>% 
  group_by(category_new) %>%
  mutate(category_new_count = paste0(category_new, '\nn = ', n())) %>% 
  ungroup()

# Simple bar plot of top variance explained
p.var.ibd <- ggplot(top.var.explained.df, aes(fill=category_new,x='IBD')) +
  geom_bar(linewidth = 0.2, colour='black',) +
  theme_classic() +
  scale_fill_manual(values = umap.category.palette,
                    name = 'Major\npopulation') +
  labs(x='', y='Number of loci')

# Simple bar plot of top variance explained
p.var.cduc <- ggplot(filter(top.var.explained.df, gwas_trait != 'IBD'),
                     aes(fill=category_new,x=gwas_trait)) +
  geom_bar(linewidth = 0.2, colour='black',) +
  theme_classic() +
  scale_fill_manual(values = umap.category.palette,
                    name = 'Major\npopulation') +
  labs(x='', y='Number of loci') 

# Build the patchwork plot
p.var.explained.bars <- (p.var.ibd | p.var.cduc) +
  plot_layout(guides='collect') &
  theme(legend.key.size = unit(0.8,"line"),
        legend.position = 'bottom',) &
  guides(fill=guide_legend(nrow=3))
p.var.explained.bars

# Count the number of annotations and make a summary barchart of those
p.var.annot <- ggplot(mutate(top.var.explained.df,
                             facet_group = ifelse(Level == 2, category_new, Level)),
                      aes(fill=category_new, y=fct_rev(fct_infreq(label_new)))) +
  geom_bar(linewidth = 0.2, colour='black',) +
  theme_classic() +
  scale_fill_manual(values = umap.category.palette,
                    name = 'Major population', guide='none') +
  scale_alpha_manual(values = c(0.5, 1), guide = 'none') +
  facet_grid(rows=vars(facet_group), scales='free', space='free') +
  labs(y='Annotations', x='Number of colocalising IBD loci') +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank())

print(p.var.annot)



# Build the patchwork plot
p.coloc.fig <- ((loci.pi.p/p.coloc.enrichment/p.var.explained.bars) | free(p.var.annot,side='b')) +
  plot_annotation(tag_levels = 'a')

if (SAVE.PLOTS == TRUE){
  ggsave(plot=p.coloc.fig,
         paste0(out.dir,"coloc_summary-PP4gt",
                gsub('[.]', 'pt', PP4.loose.threshold),".pdf"),
         width = 12, height = 10, device = cairo_pdf)
}

plot.coloc.hits.df <- known.coloc.df %>%
# plot.coloc.hits.df <- all.coloc.df %>%
  filter(PP.H4.abf > PP4.loose.threshold) %>% 
  # filter(PP.H4.abf > PP4.strict.threshold) %>%
  arrange(-PP.H4.abf) %>%  
  mutate(gene_var = paste0(symbol,'\n',qtl_lead)) %>% 
  mutate(gene_var = factor(gene_var, levels = unique(gene_var)),
         gwas_lead = factor(gwas_lead, levels = unique(gwas_lead)),
         in_OTAR_genetics = phenotype_id %in% otar.hits,
         in_OTAR_genetics_shape = case_when(phenotype_id %in% otar.hits ~ 1,
                                      TRUE ~ 10)) %>% 
  group_by(gwas_lead) %>% 
  mutate(ensgene_grouped_integer = 20+as.integer(factor(phenotype_id, levels = unique(phenotype_id)))) %>%
  ungroup()

###################
# Which colocs are transcription factors or cofactors
###################

tf.df <- read_csv(paste0(data.dir,'DatabaseExtract_v_1.01.csv')) %>%
  filter(`Is TF?` == 'Yes')
cf.df <- read_csv(paste0(data.dir,'BrowseTCOF  TcoF-DB.csv'))


# Save coloc hits for transmapping
if (SAVE.FILES == TRUE){
  # Filter for important bits
  save.coloc.genes.df <- dplyr::select(plot.coloc.hits.df, qtl_lead, condition_name, phenotype_id,
                                symbol) %>%
    distinct() %>% 
    rename(variant_id = qtl_lead)
  # Save
  write_tsv(save.coloc.genes.df,
            paste0('/home/ubuntu/eqtl/results/',tissue,'-coloc.tsv'))
  
  # Filter again for only TFs/TcoFs
  save.coloc.tf.genes.df <- filter(save.coloc.genes.df, symbol %in% c(cf.df$Symbol, tf.df$`HGNC symbol`))
  # Save
  write_tsv(save.coloc.tf.genes.df,
            paste0('/home/ubuntu/eqtl/results/',tissue,'-coloc_TFs.tsv'))
  
}

####################
# Plot colocs at a locus
####################

plot.loci.df <- known.coloc.df %>%
  group_by(annotation_type, gwas_lead, phenotype_id, gwas_trait) %>% 
  slice_max(PP.H4.abf, n=1) %>% 
  ungroup() %>% 
  arrange(condition_name, chr) %>% 
  filter(PP.H4.abf > PP4.loose.threshold)


full.trait.list <- c('AllergyAsthma'='AllergyAsthma_Zhu_NatGen2018.txt.gz',
                     'Alzheimers'='Alzheimers_Jansen_NatGen2019.txt.gz',
                     'AnkylosingSpondylitis'='AnkylosingSpondylitis_IGAS_NatGen2013.txt.gz',
                     'Asthma'='Asthma_UKBB_2018.txt.gz',
                     'AtopicExzema'='AtopicExzema_Sliz_JACI2022.txt.gz',
                     'Celiac'='Celiac_Trynka_NatGen2021.txt.gz',
                     'ChildhoodOnsetAsthma'='ChildhoodOnsetAsthma_Ferreira_AmJHGen2019.txt.gz',
                     'Coeliac'='Coeliac_Dubois_NatGen2010.txt.gz',
                     'ColonCancer'='ColonCancer_Rashkin_NatComm2020.txt.gz',
                     'CRC'='ColorectalCancer_Fernandez_NatGen2022.txt.gz',
                     'CD' = 'CrohnsDisease_DeLange_NatGen2017.txt.gz',
                     'cd' = 'cd_allchr_summary_stats_metal_nogc_Oct2023_eur_tier2_vs2_sorted.txt.gz',
                     'Diverticulosis'='Diverticulosis_UKBB_2018.txt.gz',
                     'Dupytren'='Dupytren_UKBB_2018.txt.gz',
                     'Fasciitis'='Fasciitiss_UKBB_2018.txt.gz',
                     'GravesDisease'='GravesDisease_Sakaue_NatGen2021.txt.gz',
                     'IBS'='IBS_Eijsbouts_NatGen2021.txt.gz',
                     'IBD' = 'InflammatoryBowelDisease_DeLange_NatGen2017.txt.gz',
                     'ibd' = 'ibd_allchr_summary_stats_metal_nogc_Oct2023_eur_tier2_vs2_sorted.txt.gz',
                     'MultipleSclerosis'='MS_IMSGC_NatGen2021.txt.gz',
                     'Osteoarthritis'='Osteoarthritis_Tachmazidou_NatGen2019.txt.gz',
                     'Osteoarthrosis'='Osteoarthrosis_UKBB_2018.txt.gz',
                     'Parkinsons'='Parkinsons_Nalls_LancetNeurol2019.txt.gz',
                     'PrimaryBiliaryCirrhosis'='Primary_Biliary_Cirrhosis_Liu_NatGen2012.txt.gz',
                     'Psoriasis'='Psoriasis_Tsoi_NatGen2012.txt.gz',
                     'SclerosingCholangitis'='SclerosingCholangitis_Ji_NatGen2016.txt.gz',
                     'SLE'='SystemicLupusErythematosus_Bentham_NatGen2015.txt.gz',
                     'T2D'='Type2DiabetesMellitus_Sakaue_NatGen2021.txt.gz',
                     'UC' = 'UlcerativeColitis_DeLange_NatGen2017.txt.gz',
                     'uc' = 'uc_allchr_summary_stats_metal_nogc_Oct2023_eur_tier2_vs2_sorted.txt.gz')

trait.list <- full.trait.list[unique(plot.loci.df$gwas_trait)]
trait.list <- full.trait.list[c('CD', 'UC', 'IBD')]
# trait.list <- full.trait.list[c('UC')]

# Get LD clumps
ld.min.indices <-  get_ld_clumps_per_level(clumps.path,
                                cond.ind.df, sig.interaction.sumstat.df)


plot.gene.list <- c('MAML2')

plot.gene.list <- c('MAML2', 'ZMIZ1', 'PSEN2', # Notch
                    'RASGRP1', 'LPIN3','MYC', 'FUBP1', 'FERMT1', # Wnt
                    'NDUFAF1', 'PRKCB', 'PSEN2', # Novel drug targets
                    'ITGA4', 'JAK2', 'ITGAL', 'TNFSF15', 'TNFRSF14') # Existing drug targets

# Will plot base loci in the input df that are in the plot.gene.list
for (trait in rev(names(trait.list))){
  gwas.path <- paste0('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/gwas_sumstats_final/',
                      full.trait.list[trait])
  eqtl.path <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
  plot_coloc_locus(trait,plot.loci.df,gwas.path,eqtl.path,'base','plot',
                   PP4.loose.threshold,out.dir, ld.min.indices,
                   plot.gene.list, all_variants=FALSE)
}

# Will plot all interaction loci in the input df
for (trait in rev(names(trait.list))){
  gwas.path <- paste0('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/gwas_sumstats_final/',
                      full.trait.list[trait])
  eqtl.path <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_12-multi_tissue_interaction_results/TensorQTL_eQTLS/'
  plot_coloc_locus(trait,plot.loci.df,gwas.path,eqtl.path,'interaction','plot',
                   PP4.loose.threshold,out.dir, ld.min.indices,
                   all_variants=FALSE)
}


# Make the drugs dataframe for David O
plot.coloc.hits.df %>% 
  filter(gwas_trait %in% ibd.traits) %>% 
  dplyr::select(label_new, category_new,annotation_type, Level, tissue,
                phenotype_id, gwas_lead, gwas_pval,gwas_trait, 
                qtl_lead, qtl_pval, PP.H4.abf) %>% 
  write_tsv('~/eqtl/results/coloc_results_simple.txt')


