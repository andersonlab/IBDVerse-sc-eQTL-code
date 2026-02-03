library(ggrepel)
library(shadowtext)
library(patchwork)
library(tidyverse)
library(grid)
library(ggpubr)

##############
# Variables
###############

server <- 'farm'
# server <- 'openstack'
tissue <- 'multi_tissue'


if (server == 'farm'){
  repo.dir <- paste0('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/',
                     'tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/')
  sumstats.all.basedir <- paste0('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/',
                                 'tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/')
  out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/eqtls/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  sumstats.all.basedir <- paste0('/home/rstudio/eqtl/data/IBDverse/',tissue,'/pseudobulk-base/')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/eqtls/'
  
}

data.dir <- paste0(repo.dir,'/data/')

# Plotting
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE

##################
# Read in files
##################


source(paste0(repo.dir,'qtl_plot/helper_functions.R'))


sumstat.df <- read_eqtls(sumstats.all.basedir)


sig.sumstat.df <- sumstat.df %>% 
  filter(qval < 0.05) %>% 
  mutate(annotation_type = factor(annotation_type,
                                  levels=rev(unlist(level.mapping))))


##################
# Make figure 1
##################

# Count number of egenes
df.n_egenes.sig.mean <- sig.sumstat.df %>% 
  group_by(name) %>% 
  add_count(name = 'num_eGenes') %>%
  ungroup %>% 
  distinct(annotation, tissue, Level,label_new,category_new,annotation_type,
           n_indiv, num_eGenes, mean_n_cells, mean_n_genes) %>% 
  arrange(-n_indiv)

plot.df.n_egenes.sig.mean <- df.n_egenes.sig.mean  %>% 
  filter(tissue == 'Cross-site')
  


# Add relevant plotting metadata: n_donors, average cells per donor, category
# Try alpha for 
plot_negenes_power <- function(df, x.var){
  if (x.var == 'mean_n_cells'){
    x.scale.lower <- 5
    x.scale.higher <- 7000
    x.lab <- 'Mean number of cells per donor'
    leg.pos <- c(0.8, 0.25)
  }
  if (x.var == 'mean_n_genes'){
    x.scale.lower <- 5000
    x.scale.higher <- 28000
    x.lab <- 'Mean number of genes expressed'
    leg.pos <- c(0.8, 0.25)
  }
  if (x.var == 'n_indiv'){
    x.scale.lower <- 30
    x.scale.higher <- 400
    x.lab <- 'Number of individuals'
    leg.pos <- c(0.25, 0.8)
  }
  p <- ggplot(df, aes(y = num_eGenes, x = .data[[x.var]])) +
    geom_point(aes(fill = category_new, size = annotation_type), shape = 21, stroke = 0.7) +
    scale_size_manual(values = annot.class.sizes,
                      name = 'Annotation granularity') +
    scale_x_log10(limits=c(x.scale.lower,x.scale.higher), n.breaks=6) +
    scale_y_log10(limits=c(1,20000)) +
    scale_fill_manual(values = umap.category.palette,
                      name = 'Major population') +
    theme_classic() + 
    guides(size = guide_legend(order = 1),
           fill = guide_legend(override.aes = list(size=5), order = 0,ncol=2)
    ) +
    labs(x = x.lab, y = 'Number of eGenes') +
    theme(legend.background = element_rect(fill = "transparent"),
          # legend.position = c(0.735,0.23),
          # legend.box = 'horizontal'
          # legend.position = 'none'
    )
  return(p)
}


for (var in c('mean_n_cells', 'mean_n_genes', 'n_indiv')){
  plot_negenes_power(plot.df.n_egenes.sig.mean, var)
  print(var)
  if (SAVE.PLOTS == TRUE){
    ggsave(paste0(out.dir,"n_eGenes_",var,".pdf"),
           width = 8, height = 4,device = cairo_pdf)
  }
}

# Animate
for (cur_level in 0:3){
  print(cur_level)
  cur_df <- filter(plot.df.n_egenes.sig.mean, Level < cur_level)
  print(plot_negenes_power(cur_df, 'mean_n_cells'))
  }

# Regression model

# Simple linear model but accounts for all three variables
summary(lm(log(mean_n_cells) ~ log(n_indiv) + log(mean_n_genes) + log(mean_n_cells) ,
   data = filter(plot.df.n_egenes.sig.mean,
                 Level == 2)))

##################
# Make figure 2
##################

df.count.sig.egenes <- sig.sumstat.df %>% 
  mutate(Gene_class = factor(case_when(biotype == 'protein_coding' ~ 'Protein coding',
                                       biotype == 'lncRNA' ~ 'lncRNA',
                                       TRUE ~ 'Other'),
                             levels = c('Other', 'lncRNA', 'Protein coding'))) %>% 
  arrange(Level) %>% 
  distinct(phenotype_id, .keep_all = TRUE) 

# Make df to keep track of tissue exclusivity per eGene
df.count.sig.egenes.tissue <- sig.sumstat.df %>% 
  mutate(tissue = factor(tissue, 
                         levels = rev(unname(tissue.mapping)))) %>%
  group_by(phenotype_id) %>%
  summarise(Tissue_detection = paste(sort(unique(tissue)), collapse = ", "),
            .groups = "drop" ) %>% 
  mutate(Tissue_detection = str_remove(Tissue_detection, 
                                       paste(',', unname(tissue.mapping[1]))),
         Tissue_detection_count = str_count(Tissue_detection, ','),
         Tissue_detection_count_label = factor(case_when(
           Tissue_detection == unname(tissue.mapping)[1] ~ ' ',
           Tissue_detection_count == 0 ~ 'One site',
           Tissue_detection_count == 1 ~ 'Two sites',
           TRUE ~ 'All sites'),
           levels = c(' ',
                      'All sites', 
                      'Two sites', 
                      'One site')))

# Bar plot of tissue sharing of eGenes
ggplot(df.count.sig.egenes.tissue,aes(y=Tissue_detection, fill=Tissue_detection)) +
  geom_bar(linewidth = 0.2, colour='black') +
  scale_fill_manual(values=tissue.palette, guide='none') +
  theme_classic() +
  xlim(0,6500) +
  labs(y = '', x = 'Number of eGenes') +
  facet_grid(rows=vars(Tissue_detection_count_label),scales = 'free', space = 'free',
             switch ='y') +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"n_eGenes_tissue_sharing.pdf"),
         width = 4.5, height = 4,device = cairo_pdf)
}


# How many eGenes unique to an annotation do we find?
in.n.annots.breaks <- list(  
  breaks = c(0, 1, 2, 100),
  labels = c("1", "2", "3 or more")
)
in.n.annots <- sumstat.df %>% 
  filter(qval < 0.05) %>%
  # filter(bonf_pvalue < 0.05) %>%
  distinct(phenotype_id, label_new,.keep_all = TRUE) %>% 
  group_by(phenotype_id) %>% 
  mutate(egene_in_n_annotations = n(),
         everything_present = any(label_new == level.mapping[['0']])) %>%
  ungroup() %>% 
  dplyr::select(phenotype_id, symbol, egene_in_n_annotations,
                label_new, everything_present, category_new, Level) %>% 
  arrange(desc(egene_in_n_annotations)) %>% 
  distinct(phenotype_id,.keep_all = TRUE) %>% 
  mutate(binned_egene_in_n_annotations = cut(
    egene_in_n_annotations,
    breaks = in.n.annots.breaks$breaks,
    labels = in.n.annots.breaks$labels,
    include.lowest = TRUE
  ))

# Plot how many eGenes total are captured and how many the Unannotated category captures alone
plot_n_egenes_count <- function(plot_df, fill_col){
  if (fill_col == 'annotation_type'){
    palette <- annot.class.palette
    y.axis.title <- 'Granularity'
    legend.title <- 'Lowest granularity eGene\ncan be observed at'
  }
  if (fill_col == 'Gene_class'){
    palette <- c('lncRNA'='#fc8d62', 'Protein coding'='#e78ac3','Other' = "grey")
    y.axis.title <- 'Type'
    legend.title <- 'Transcript class'
  }
  if (fill_col == 'Tissue_detection'){
    palette <- tissue.palette
    y.axis.title <- 'Tissue'
    legend.title <- ''
  }
  if (fill_col == 'binned_egene_in_n_annotations'){
    palette <- 'Reds'
    y.axis.title <- 'Sharing'
    legend.title <- 'Number of annotations\neGene is observed in'
  }
  p <- ggplot(plot_df,
         aes(x='', fill=.data[[fill_col]])) +
    geom_bar(linewidth = 0.25, colour='black') +
    # geom_text(stat = 'count',aes(label = after_stat(count),
    #                              y=after_stat(count)),colour='black',
    #           position = position_stack(vjust = 0.5),
    #           size=4.8,fontface='bold',) +
    theme_classic() +
    coord_flip() +
    # scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ylim(0,21000) +
    labs(x = y.axis.title, y = 'Number of eGenes') +
    theme(axis.ticks.y=element_blank()) +
    guides(fill=guide_legend(nrow=3,reverse = TRUE))
  if (fill_col != 'Tissue_detection'){
  }
  if (fill_col != 'binned_egene_in_n_annotations'){
    p <- p + scale_fill_manual(values = palette, name = legend.title) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x = element_blank())
  }
  if (fill_col == 'binned_egene_in_n_annotations'){
    p <- p +  scale_fill_brewer(palette = palette, name=legend.title, direction=-1)
  }
  return(p)
}

p2.top <- plot_n_egenes_count(df.count.sig.egenes, 'annotation_type')
p2.middle <- plot_n_egenes_count(df.count.sig.egenes, 'Gene_class')
p2.bottom <- plot_n_egenes_count(left_join(df.count.sig.egenes,
                                           in.n.annots, by='phenotype_id'),
                                 'binned_egene_in_n_annotations')
# Animate
for (cur_level in 1:3){
  print(cur_level)
  cur_df <- filter(df.count.sig.egenes, Level < cur_level)
  print(plot_n_egenes_count(cur_df, 'annotation_type') + theme_classic())
}

p2.top/p2.middle/p2.bottom + plot_layout(guides = 'collect')

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"n_eGenes_annotation.pdf"),
         width = 14, height = 4.5,device = cairo_pdf)
}

# Plot how many annots each eGene is found in
ggplot(in.n.annots, aes(x=egene_in_n_annotations, fill=everything_present)) + 
  geom_bar() + 
  theme_classic() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = c(0,max(in.n.annots$egene_in_n_annotations))) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = "black")) +
  # scale_fill_manual(values=c('TRUE'='grey25',
  #                            'FALSE'='grey65')) +
  labs(x='Number of annotations an eGene is found in',
       y='Number of eGenes') +
  guides(fill=guide_legend(title="eGene is found in\nAll Cells"))

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"egene-found_in_n_annots.pdf"),
         width = 7, height = 5.5,device = cairo_pdf)
}

in.one.annot <- in.n.annots %>% 
  filter(egene_in_n_annotations == 1) %>% 
  group_by(label_new) %>% 
  mutate(total_egenes = n()) %>% 
  ungroup() 

# Count the number of times an exclusive eGene is found per annotation
ggplot(filter(in.one.annot, Level != 0, total_egenes > 5),
       aes(y=forcats::fct_infreq(label_new),fill=category_new)) +
  geom_bar() +
  scale_fill_manual(values=umap.category.palette, guide='none') +
  theme_classic() +
  # theme(panel.grid.major.y = element_blank(),
  #       panel.grid.minor.y = element_blank()) +
  labs(y = '', x = 'Number of eGenes unique to annotation') +
  # xlim(c(0,400))+
  facet_grid(rows=vars(-Level),scales = 'free', space = 'free') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"egene-found_in_1_annots.pdf"),
         width = 6.5, height = 6.5,device = cairo_pdf)
}


# Venn diagram of the tissue sharing of eGenes
tissue.frequency <- df.count.sig.egenes.tissue %>%
  count(Tissue_detection) %>% 
  deframe() 

# What proportion are only found in a single site
(tissue.frequency[['Blood']] + tissue.frequency[['Rectum']] + tissue.frequency[['Terminal ileum']])/sum(tissue.frequency)

# How many are only found in cross-site
(tissue.frequency[['Cross-site']])/sum(tissue.frequency)

##################
# Make figure 4
##################
# Effect sizes of lead eQTLs for each annotation and location relative to TSS


# Distribution of effect size stratified by granularity
f4.p1 <- ggplot(filter(sig.sumstat.df, phenotype_id %in% in.one.annot$phenotype_id),
                       aes(x=abs(slope), fill=annotation_type,
                           colour=annotation_type)) + 
  geom_density(alpha=0.7,) +
  # geom_boxplot(alpha = 0.4) + 
  theme_classic() + 
  scale_fill_manual(values = annot.class.palette) +
  scale_colour_manual(values = annot.class.palette) +
  labs(x='Abs(β)') +
  guides(fill=guide_legend(title="Granularity"),
         color='none') + 
  theme(legend.position = c(0.8, 0.7))

maf.sig.sumstat.df <- sig.sumstat.df %>% 
  mutate(maf = case_when(af > 0.5 ~ 1-af,
                         TRUE ~ af))

# Correlating MAF to effect size
f4.p2 <- ggplot(maf.sig.sumstat.df, aes(x=maf, y=abs(slope),colour=annotation_type,
                               group=annotation_type)) +
  geom_point(alpha=0.15)  + 
  geom_smooth(method = lm, size = 3, color = "black",,se=FALSE) + 
  geom_smooth(method=lm,size=2,,se=FALSE) +
  theme_classic() +
  scale_colour_manual(values = annot.class.palette) +
  labs(y='Abs(β)', x = 'Minor allele frequency') +
  guides(colour=guide_legend(title="Granularity")) + 
  theme(legend.position = c(0.8, 0.7))

f4.p2 #+ plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"eqtl-lead_effect_granularity.pdf"),
         width = 4, height = 4, device = cairo_pdf)
}

# Does the average effect size differ between annotation (only check one resolution level)
f4.p3 <- 
  # Using median
  filter(sig.sumstat.df, Level == 1) %>%
  mutate(category_new = fct_reorder(category_new, slope, .fun='median')) %>%
  ggplot(aes(x=reorder(category_new, abs(slope)),
             y=abs(slope), colour=category_new)) + 
  geom_jitter(alpha=0.25) +
  geom_violin(alpha=0.8,draw_quantiles = TRUE)+
  geom_boxplot(width=.2, alpha=0.7,outlier.shape = NA) + 
  scale_colour_manual(values=umap.category.palette) +
  theme_classic() +
  labs(y='Abs(β)',x='') +
  guides(colour='none')

f4.p3

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"eqtl-lead_effect_granularity_by_category.pdf"),
         width = 7, height = 5.5,device = cairo_pdf)
}

##################
# Make figure 5
##################

cond.ind.df <- read_conditional_eqtls(sumstats.all.basedir)

plot.cond.ind.df <- cond.ind.df %>%
  filter(rank != 1) %>% 
  mutate(rank = case_when(rank > 2 ~ 'Tertiary or greater',
                          TRUE ~ 'Secondary')) %>% 
  group_by(annotation_type,category_new, label_new, rank) %>%
  mutate(n_eQTLs=sum(abs(direction))) %>%
  ungroup() %>%
  arrange(-n_eQTLs) %>%
  mutate(label_new = factor(label_new, levels = unique(label_new)),
         category_new = factor(category_new, levels = unique(category_new)),
         annotation_type = factor(annotation_type, levels = unique(annotation_type))) %>% 
  distinct(label_new,n_eQTLs,rank, .keep_all = TRUE)


# Plot number of conditionally independent SNPs
f5.p1 <- ggplot(plot.cond.ind.df,
       aes(x=label_new, y=n_eQTLs, alpha=rank, fill=category_new)) +
  geom_bar(stat = 'identity') +
  facet_grid(~annotation_type, scales = 'free', space = 'free') +
  geom_hline(yintercept = 0) +
  geom_text(data = filter(plot.cond.ind.df, rank == 'Secondary'),
            aes(label = abs(n_eQTLs)), colour='black', angle=90, hjust=0,
            size=3, vjust=1,show.legend = FALSE) +
  geom_text(data = filter(plot.cond.ind.df, rank != 'Secondary'),
            aes(label = abs(n_eQTLs)), colour='black', angle=90, hjust=0,
            size=3, vjust=0,show.legend = FALSE) +
  scale_alpha_discrete(range = c(0.5, 1),
                       name = 'eQTL rank') +
  theme_classic() +
  scale_fill_manual(values = umap.category.palette,
                    name = 'Major category\ncolour',
                    guide='none') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = c(0.9, 0.9),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        # panel.grid.major.y = element_line(color = "lightgrey",
        #                                   size = 0.25,
        #                                   linetype = 1),
        panel.grid.minor.x = element_blank())+
  labs(x = '', y = 'Number of non-lead\nindependent eQTLs',
       title= '')
f5.p1

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"conditional-eQTLs.pdf"),
         width = 15, height = 6,device = cairo_pdf)
}

f5.p2 <- ggplot(filter(cond.ind.df,Level != 2),
# f5.p2 <- ggplot(cond.ind.df,
                aes(y=abs(start_distance), fill=as.factor(rank),x=as.factor(rank))) +
  geom_jitter(aes(colour=as.factor(rank)),alpha=0.25) +
  geom_boxplot(outlier.shape = NA, alpha=0.5, colour='black') +
  # stat_compare_means(label.x=2,label.y=900000,size=3,
                     # method='kruskal.test',label = "p.format") +
  stat_compare_means(comparisons = list(c(1,2)), label.y = c(7e5),size=3,label = "p.signif") +
  stat_compare_means(comparisons = list(c(1,3)), label.y = c(8.5e5),size=3,label = "p.signif") +
  facet_grid(~category_new,space = 'free_y') +
  theme_classic() +
  scale_fill_viridis_d(direction = -1,
                       guide = 'none',
                       name = 'Rank') +
  scale_colour_viridis_d(direction = -1,
                         guide = 'none') + 
  theme(strip.background = element_blank()) +
  labs(x = 'eQTL rank', y = 'Distance from TSS (bp)')

f5.p1 / f5.p2 #+ plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

lm(data=filter(cond.ind.df,annotation_type != 'category__machine'),formula = 'abs(start_distance) ~ rank + category_new')

# Plot the multiple eQTL distances from TSS for a particular gene
cur_gene <- 'MAML2'
ggplot(filter(cond.ind.df, symbol == cur_gene)) +
  geom_point(aes(x=position, fill=tissue,y=label_new),shape=21, size=5, alpha=0.75) +
  geom_rect(data=grch38 %>% filter(symbol == 'MAML2'),aes(ymin=0,ymax=-1,
                                                          xmin=start, xmax=end),
            fill='black') +
  scale_fill_manual(values = tissue.palette) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  labs(y = '', x = 'Genomic position (bp)')


#############
# Check vs TFs
##############

human.tf.df <- read_csv('/nfs/users/nfs_o/oa3/eqtl/data/TF_list_lambert_2017.csv')

#############
# Supp fig
##############

# heatmap of eGenes df
df.heatmap.sig.egenes <- sig.sumstat.df %>%
  filter(qval < 0.05)

# Heatmap of eGenes
ggplot(df.heatmap.sig.egenes, aes(x=phenotype_id, y=label_new, fill=category_new)) + 
  geom_tile() +
  theme_bw()

