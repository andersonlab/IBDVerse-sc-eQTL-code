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
  repo.dir <- '~/eqtl/code/IBDVerse-sc-eQTL-code/'
  interaction.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_12-multi_tissue_interaction_results/TensorQTL_eQTLS/'
  out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/interaction/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code'
  interaction.basedir <- paste0('/home/rstudio/eqtl/data/IBDverse/',tissue,'/pseudobulk-interaction/')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/interactions/'
}
data.dir <- paste0(repo.dir,'/data/')

# Plotting
SAVE.PLOTS <- TRUE

##############
# Functions
##############

source(paste0(repo.dir,'qtl_plot/helper_functions.R'))

#############
# Code
#############

interaction.sumstat.df <- read_ieqtls(interaction.basedir)


sig.interaction.sumstat.df <- interaction.sumstat.df %>% 
  filter(pval_adj_bh < 0.05,
         tissue != 'Cross-site',
         interaction != 'ses_inflamed')

###########
# Make figure 1
############

# Make a bar chart of number of interactions

f1.p1 <- ggplot(sig.interaction.sumstat.df, aes(x = category_new,
                                    fill = category_new)) + 
          geom_bar() + 
          theme_classic() +
          labs(y = 'Number of ieGenes', x = '') +
  scale_fill_manual(values = umap.category.palette,
                    name = 'Major\nPopulation') + 
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(interaction_new~tissue, scales = 'free_y') 
f1.p1

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"n_ieGenes.pdf"),
         width = 5.5, height = 5.5,device = cairo_pdf)
}


############
# Figure 3
############

# Plot top 3 (by nES) terms from medication, inflammation (mod/sev), sex and age
positive_ES_terms <- sig.gsea.df %>% 
  filter(interaction == 'ses_cd_binned',
         category_new %in% c('Enterocyte', 'Stem')) %>%
  mutate(annot_id_clean = str_to_sentence(gsub('_|REACTOME', ' ', annot_id))) %>% 
  filter(annot_coef > 0)

ggplot(positive_ES_terms, aes(y=annot_id_clean, size = annot_coef, x=-log10(qvalue_bh),
                         colour=category_new,label=interaction)) +
  geom_point(alpha=0.75) + 
  # geom_text_repel(size=3) +
  scale_colour_manual(values = umap.category.palette,guide='none')  +
  theme_bw() +
  labs(size = 'Normalised enrichment score', y='Geneset',
       x=expression(-log[10](q[FGSEA])), colour='Major category')+
  # scale_size(breaks=c(0,1,2,3,4)) +
  theme(panel.grid.major.y = element_blank(),
        legend.text = element_text(size=12),
        legend.title.align = 0,
        legend.position = "bottom", 
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.title = element_text(size=14, vjust = .5, hjust = .3))

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,"ieGene_GSEA.pdf"),
         width = 9.6, height = 3.5,device = cairo_pdf)
}

# Which leading edge genes are significant in the ieGenes
leading_edge_genes <- positive_ES_terms$leadingEdge %>% 
  strsplit(split = ',') %>% 
  unlist() %>% 
  unique() %>% 
  sort()

# Filter sig.sumstat.df for leading edge genes
leading.edge.sig.sumstat.df <- sig.interaction.sumstat.df %>% 
  filter(symbol %in% leading_edge_genes)

# Filter sig.gsea.df for rows where at least one significant ieGene is in the leading edge
sig.iegene.sig.gsea.df <- sig.gsea.df %>%
  filter(sapply(strsplit(as.character(leadingEdge), ","), 
                function(x) any(trimws(x) %in% sig.interaction.sumstat.df$symbol)))

############
# Figure 4
############
plot_gsea_volcana <- function(cur_gsea_row,interaction.sumstat.df, only.sig = FALSE){
    print(cur_gsea_row$label_new)
    print(cur_gsea_row$interaction_new)
    plot.df <- interaction.sumstat.df %>%
      filter(annotation == cur_gsea_row$annotation,
             interaction == cur_gsea_row$interaction)

    
    leadingedge <- strsplit(cur_gsea_row$leadingEdge, split=',')[[1]]
    print(leadingedge)
    if (only.sig){
      only.sig.df <- filter(plot.df, symbol %in% filter(sig.interaction.sumstat.df, 
                                                    label_new == cur_gsea_row$label_new,
                                                    symbol %in% leadingedge,
                                                    interaction_new == cur_gsea_row$interaction_new)$symbol)
      if (nrow(only.sig.df) == 0){
        return("None")
      }
    }
    
    plot.text <- paste0('\n\nFGSEA q-value: ',
                        round(cur_gsea_row$qvalue_bh,5), '\n',
                        'FGSEA normalised enrichment score: ',
                        round(cur_gsea_row$annot_coef,2))
    
    p <- ggplot(filter(plot.df, symbol %in% leadingedge),
                aes(x=b_gi, y=-log10(pval_gi),
                    fill=category_new)) +
                  # colour=symbol %in% leadingedge)) +
      geom_point(data=filter(plot.df, !symbol %in% leadingedge), alpha=0.2,
                 shape = 21, colour='white') +
      geom_point(size=3, shape = 21, colour='black') +
      geom_point(data=filter(plot.df, symbol %in% filter(sig.interaction.sumstat.df, 
                                                         label_new ==cur_gsea_row$label_new,
                                                         symbol %in% leadingedge,
                                                         interaction_new == cur_gsea_row$interaction_new)$symbol),
                 size=1.5, shape = 21, fill='black') +
      annotate("text",  x=Inf, y = Inf, label = plot.text, vjust=1, hjust=1,
               size=3) +
      geom_label_repel(aes(label=symbol), size=3, max.overlaps = 30, fill='white') +
      theme_classic() +
      scale_fill_manual(values = umap.category.palette) +
      theme(legend.position = 'none',
            axis.title = element_text(size=7),
            axis.text = element_text(size=5),) +
      labs(x = 'Genotype:Interaction β',
           y=expression(-log[10](p[G:I])),
           title=paste0(cur_gsea_row$label_new,'\n', cur_gsea_row$interaction_new, 
                       '\n',cur_gsea_row$annot_id_clean)) 
    return(p)
}

plot_top_terms <- positive_ES_terms %>% 
  arrange(annot_coef) %>% 
  filter(interaction_new == 'Inflammation') %>% 
  distinct(leadingEdge, .keep_all = TRUE) %>% 
  distinct(annot_id, .keep_all = TRUE)

volcano_plot_list <- list()

for (i in 1:nrow(plot_top_terms)){
  volcano_plot_list[[i]] <- plot_gsea_volcana(plot_top_terms[i,], interaction.sumstat.df, only.sig = TRUE)
  print(i)
  print(volcano_plot_list[[i]])
}
# patchwork::wrap_plots(volcano_plot_list, nrow = 2)

# Just metal related terms
metal_terms <- positive_ES_terms %>%  
  filter(str_detect(annot_id_clean, 'etal'))

volcano_metal_plot_list <- list()
for (i in 1:nrow(metal_terms)){
  print(plot_gsea_volcana(metal_terms[i,], interaction.sumstat.df))
}
# patchwork::wrap_plots(volcano_metal_plot_list, nrow = 2)

# Just terms where at least one of the leading edge genes is in my sig.sumstat.df
plot.sig.iegene.sig.pos.terms <- positive_ES_terms %>% 
  filter(sapply(strsplit(as.character(leadingEdge), ","), 
                function(x) any(trimws(x) %in% sig.interaction.sumstat.df$symbol))) %>% 
  arrange(annot_coef) %>%
  distinct(leadingEdge, .keep_all = TRUE) %>% 
  distinct(annot_id, .keep_all = TRUE)

for (i in 1:nrow(plot.sig.iegene.sig.pos.terms)){
  print(plot_gsea_volcana(plot.sig.iegene.sig.pos.terms[i,], interaction.sumstat.df, only.sig = TRUE))
}

############
# Interaction colocs
############

interaction.colocs <- get_colocs('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed/')

sig.colocs <- interaction.colocs %>% 
  filter(PP.H4.abf > 0.75)

# Filter sig.colocs for leading edge genes
leading.edge.sig.colocs <- sig.colocs %>% 
  filter(symbol %in% leading_edge_genes)

ggplot(sig.colocs, aes(fill=category_new, x=gwas_lead_pos, y=PP.H4.abf, size=-log10(gwas_pval))) + 
  geom_point(shape=21) + 
  geom_text_repel(aes(label=gwas_trait), size=3, max.overlaps = 30, fill='white') +
  theme_classic() + 
  scale_x_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = umap.category.palette) +
  facet_grid(~chr, scales = 'free', space='free') +
  labs(x = 'GWAS lead position', y = 'Posterior probability of colocalisation')
  
  
  
  
  