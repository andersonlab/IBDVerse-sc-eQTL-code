library(plyr)
library(RColorBrewer)
library(tidyverse)

##############
# Variables
###############

server <- 'farm'
tissue <- 'multi_tissue'
cohort <- 'atlas'
# cohort <- 'eqtl'

if (server == 'farm'){
  repo.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/'
  metadata.basedir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2024_12_31-multi_tissue_base_results/metadata/'
  out.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  metadata.basedir <- paste0('/home/rstudio/eqtl/data/IBDverse/',tissue,'/metadata/')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/atlas/'
  
}

data.dir <- paste0(repo.dir,'/data/')

# Plottings
SAVE.PLOTS <- TRUE
SAVE.FILES <- FALSE

##############
# Functions
##############

source(paste0(repo.dir,'/qtl_plot/helper_functions.R'))

##############
# Main
###############

# Not filtering for cells<5 per sample-annot combination
if (cohort == 'eqtl'){
cell.metadata <- read_csv(paste0(metadata.basedir,'/',cohort,'_processed.obs.csv')) %>% 
  dplyr::select(-tissue) %>% 
  left_join(annot.mapping, by=c('predicted_labels_tissue' = 'label_machine'))
} else{
  cell.metadata <- read_csv(paste0(metadata.basedir,'/',cohort,'_processed.obs.csv')) %>% 
    mutate(label_machine = paste0(leiden,'_', tissue)) %>% 
    dplyr::select(-tissue) %>% 
      left_join(annot.mapping, by='label_machine') %>% 
    drop_na(label_new)
}



plot_cell_counts <- function(df, position_var, fill_var,legend_pos='left',
                             x_lab='',y_max=max(table(cell.metadata$label_new))){
  guide_nrow <- 3
  if (fill_var == 'tissue'){
    palette <- tissue.palette
    legend.title <- 'Anatomical site'
  }
  if (fill_var == 'disease_status'){
    palette <- disease.palette
    legend.title <- 'Disease status'
    guide_nrow <- 2
  }
  if (fill_var == 'inflammation_status'){
    
    legend.title <- 'Inflammation status'
    guide_nrow <- 2
  }
  if (fill_var == 'category_new'){
    palette <- umap.category.palette
    legend.title <- 'Major population'
  }
  if (fill_var == 'sex'){
    palette <- c('M'='#bebada', 'F' = '#8dd3c7')
    legend.title <- 'Sex'
    guide_nrow <- 2
  }
  if (fill_var == 'age'){
    df['age_binned'] <- round_any(df$age, 10) 
    age_labels <- c('20'='18-24','30'='25-34','40'='35-45','50'='45-54','60'='55-64','70'='65+')
    df <- df %>% mutate(age_binned_labels = age_labels[as.character(round_any(df$age, 10))])
    df['age_binned_labels'] <- factor(df$age_binned_labels, 
                                      levels = unname(age_labels))
    palette <- brewer.pal(length(age_labels),"RdYlBu")
    names(palette) <- unname(age_labels)
    fill_var <- 'age_binned_labels'
    legend.title <- 'Age'
  }
  p <- ggplot(df, aes(fill=.data[[fill_var]],x=forcats::fct_infreq(label_new))) + 
    geom_bar(position=position_var) +
    theme_classic() +
    facet_grid(~category_new, scales = 'free', space='free') +
    labs(x=x_lab,y='Proportion of cells') +
    ylim(0,y_max) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position = legend_pos,
          strip.background = element_blank()) +
    scale_fill_manual(values = palette, name=legend.title,) + 
    guides(fill=guide_legend(nrow=guide_nrow)) +
    theme(strip.text = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if (position_var == 'stack'){
    p <- p + ylab('Number of cells')
  } 
  if (position_var == 'fill'){
    p <- p + scale_y_continuous(breaks = c(0,0.5,1))
  }
  
  return(p)
}

(plot_cell_counts(cell.metadata, 'stack', 'category_new', legend_pos='none') +
    theme(strip.text = element_text(face = 'bold'))) /(
      (plot_cell_counts(cell.metadata, 'fill', 'tissue')) / 
        (plot_cell_counts(cell.metadata, 'fill', 'disease_status')) /
        # (plot_cell_counts(cell.metadata, 'fill', 'sex')) /
        (plot_cell_counts(cell.metadata, 'fill', 'age') + theme(axis.ticks.x = element_line(),
                                                                axis.text.x = element_text(angle=45,
                                                                                           hjust=1,
                                                                                           vjust=1))) +
        plot_layout(guides='collect',axis_titles = 'collect') & theme(legend.position = 'left')) +
  plot_layout(heights = c(1,1.5))

if (SAVE.PLOTS == TRUE){
  ggsave(paste0(out.dir,cohort,"-celltype_distribution.pdf"),
         width = 18, height = 8,device = cairo_pdf)
}

# Plot on a per cell type level
for (cur_annot in unique(annot.mapping$label_new)){
  print(cur_annot)
  cur.cell.metadata <- filter(cell.metadata, label_new == cur_annot) %>% 
    mutate(disease_status = factor(disease_status,
                                   levels = c('Healthy','CD')))
  if (nrow(cur.cell.metadata) == 0){
    next
  }
  plot_cell_counts(cur.cell.metadata, 'stack', 'tissue',x_lab = 'Cross-site',y_max=nrow(cur.cell.metadata)) +
    (plot_cell_counts(filter(cur.cell.metadata, tissue == 'Terminal ileum'),
                      'stack', 'disease_status',x_lab = 'TI',y_max=nrow(cur.cell.metadata)) + 
       theme(axis.ticks.y = element_blank(),
             axis.text.y = element_blank(),
             axis.title.y = element_blank())) +
    plot_annotation(title = cur_annot) +
    plot_layout(guides='collect') & 
    guides(fill = guide_legend(title.position = "top",ncol=1)) 
  
  if (SAVE.PLOTS == TRUE){
    cur_annot <- str_replace_all(cur_annot, ' ', '_')
    cur_annot <- str_replace_all(cur_annot, '/', '_')
    ggsave(paste0(out.dir,"within_annot/celltype_distribution_",cur_annot,".pdf"),
           width = 3.8, height = 2.8,device = cairo_pdf)
  }
}

# Rank cell types by a combination of how disease and tissue-specific they are
celltype.frequencies <- cell.metadata %>% 
  group_by(tissue, label_new) %>% 
  add_count(name = 'num_tissue_celltype') %>%
  ungroup() %>% 
  group_by(disease_status, label_new) %>% 
  add_count(name = 'num_disease_celltype') %>%
  ungroup() %>% 
  group_by(label_new) %>% 
  dplyr::mutate(freq_tissue_celltype = round(num_tissue_celltype/n(),3),
                freq_disease_celltype = round(num_disease_celltype/n(),3)) %>% 
  ungroup() %>% 
  dplyr::select(label_new,tissue,freq_tissue_celltype,disease_status,freq_disease_celltype) %>% 
  distinct() %>% 
  filter(freq_disease_celltype > 0.25, freq_tissue_celltype > 0.25) %>% 
  arrange(desc(tissue), disease_status, desc(freq_tissue_celltype + freq_disease_celltype))


celltype.ranking <- celltype.frequencies %>% 
  distinct(label_new,.keep_all = TRUE) %>% 
  mutate(rank = row_number())
