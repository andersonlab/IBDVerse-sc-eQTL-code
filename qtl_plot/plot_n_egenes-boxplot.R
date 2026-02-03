library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggpubr)

##############
# Variables
###############

server <- 'farm'
# server <- 'openstack'
tissue <- 'multi_tissue'


if (server == 'farm'){
  repo.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/'
  sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
  coloc.base.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed'
  coloc.interaction.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed'
  out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/'
  vcf.path <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/imputed.vcf.gz'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  sumstats.all.basedir <- paste0('/home/rstudio/eqtl/data/IBDverse/',tissue,'/pseudobulk-base/')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/coloc/'
  
}

data.dir <- paste0(repo.dir,'/data/')
normalised <- TRUE

# Plotting
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE

##################
# Read in files
##################


source(paste0(repo.dir,'qtl_plot/helper_functions.R'))

if (!exists('sumstat.df')){
  sumstat.df <- read_eqtls(sumstats.all.basedir)
}

 sig.sumstat.df <- sumstat.df %>% 
  filter(qval < 0.05) %>% 
  mutate(annotation_type = factor(annotation_type,
                                  levels=rev(unlist(level.mapping))))

if (!exists('known.coloc.df')){
  print('Loading coloc results')
  known.base.coloc.df <- get_colocs(coloc.base.dir, known_ibd_only = TRUE)
  known.interaction.coloc.df <- get_colocs(coloc.interaction.dir, known_ibd_only = TRUE) %>% 
    mutate(interaction = str_split_fixed(condition_name, "__", 3)[,3])
  
  known.coloc.df <- bind_rows(mutate(known.base.coloc.df, interaction = 'base'),
                              known.interaction.coloc.df)
}
 
# Load vcf
if (!exists('vcf.df')){
  print('Loading VCF')
  vcf.df <- read.table(vcf.path,skip=48, header=TRUE, comment.char = "")
  }


int.expression.paths <- list.files(path = paste0(sumstats.all.basedir, "../norm_data/"), 
                                 pattern = "normalised_phenotype.tsv$",
                                 recursive = TRUE,  full.names = TRUE)

raw.expression.paths <- list.files(path = paste0(sumstats.all.basedir, "../aggregated_counts/"), 
                                 pattern = "phenotype_file.tsv$",
                                 recursive = TRUE,  full.names = TRUE)

# Must have the colocs loaded first
top.coloc.snps.df <- filter(known.coloc.df, PP.H4.abf > 0.75, gwas_trait %in% ibd.traits) %>%
                            group_by(loci_ID,phenotype_id) %>% slice_max(PP.H4.abf)# %>% 
  # mutate(snp_id = str_replace_all(snp_id, '_', ':'))
top.coloc.snps <- unique(top.coloc.snps.df$snp_id)

# Pre-split known.coloc.df and vcf.df for fast lookups by variant
filtered.known.coloc.df <- filter(known.coloc.df, 
                                  snp_id %in% top.coloc.snps,
                                  gwas_trait %in% ibd.traits) %>% 
  mutate(snp_id = str_replace_all(snp_id, '_', ':')) %>% arrange(-PP.H4.abf)
coloc_lookup <- split(filtered.known.coloc.df, filtered.known.coloc.df$snp_id)

replaced_top_coloc_snps <- top.coloc.snps %>% 
  str_replace_all('_', ':')

filter.vcf.df <- filter(vcf.df, ID %in% replaced_top_coloc_snps)
vcf_lookup <- split(filter.vcf.df, filter.vcf.df$ID)


# Use a list to store results, then bind at the end
var_explained_summary_list <- list()
result_index <- 1

for (cur.annot in unique(annot.mapping$label_machine)) {
  print(cur.annot)
  # Skip if we've already processed this annotation
  if (cur.annot %in% map_chr(var_explained_summary_list, "annotation")) {
    next
  } else if (!any(grepl(cur.annot, int.expression.paths))){
    print(paste('Skipping', cur.annot))
    next
  }
  

  cur.base.annot.exp.df <- read.table(str_subset(raw.expression.paths, cur.annot),header = TRUE) %>% 
    dplyr::rename(ensgene = ENS)
  
  cur.int.annot.exp.df <- read.table(str_subset(int.expression.paths, cur.annot), header = TRUE, row.names = NULL) %>% 
    dplyr::rename(ensgene = row.names)
  
  for (cur.var in replaced_top_coloc_snps) {
    # Look up the sig.sumstat entry; skip if not found
    if (!cur.var %in% names(coloc_lookup)) next
    cur_eqtl_of_interest <- coloc_lookup[[cur.var]][1, ]
    
    for (cur.gene.name in unique(cur_eqtl_of_interest$phenotype_id)) {
      # Look up and process VCF data; skip if not found
      if (!cur.var %in% names(vcf_lookup)) next
      cur.variant.df <- vcf_lookup[[cur.var]] %>%
        pivot_longer(-c(1:9), names_to = 'genotype_id') %>%
        mutate(genotype_id = str_remove(genotype_id, '^X')) %>%
        separate(value, into = c('GT','DS','HDS','GP'), sep = ':') %>%
        mutate(Genotype = str_count(GT, "1")) %>%
        select(ID, genotype_id, Genotype, GT) %>%
        drop_na()
      # Check if the gene is in the expression data
      if (!cur.gene.name %in% cur.int.annot.exp.df$ensgene) {
        next
      }
      
      # Process gene expression for the gene of interest
      cur.gene.df <- cur.int.annot.exp.df %>%
        filter(ensgene == cur.gene.name) %>%
        pivot_longer(-ensgene, names_to = 'annotation_individual') %>%
        separate(annotation_individual,
                into = c('annotation_type','annotation','individual'),
                sep = '[.]') %>%
        mutate(
          part2 = str_split(individual, "_", simplify = TRUE)[,2],
          part3 = str_split(individual, "_", simplify = TRUE)[,3],
          genotype_id = if_else(part2 == part3, part2, paste(part2, part3, sep = "_"))
        ) %>% 
        select(-part2, -part3, -individual) %>%
        left_join(cur.variant.df, by = 'genotype_id') %>%
        drop_na()
      # Fit a linear model if there are rows and number of unique genotoypes > 1
      if (nrow(cur.gene.df) > 0 & length(unique(cur.gene.df$Genotype)) > 1) {
        lm.summary <- summary(lm(value ~ Genotype, data = cur.gene.df))
        cur.rsq <- lm.summary$adj.r.squared
        cur.beta <-  lm.summary$coefficients["Genotype", "Estimate"]
        cur.beta.error <-  lm.summary$coefficients["Genotype", "Std. Error"]
        cur.pvalue <-  lm.summary$coefficients["Genotype", "Pr(>|t|)"]
      } else {
        cur.rsq <- NA      
        cur.beta <- NA
        cur.beta.error <- NA
        cur.pvalue <- NA
      }
      
      # Extract pseudobulk expression for the gene, compute its mean and variance
      cur.base.annot.exp.gene.df <- cur.base.annot.exp.df %>%
        filter(ensgene == cur.gene.name) %>%
        select(-ensgene) %>%
        t() %>% 
        .[, 1]
      
      cur.mean.exp <- mean(cur.base.annot.exp.gene.df)
      cur.median.exp <- median(cur.base.annot.exp.gene.df)
      cur.var.exp <- var(cur.base.annot.exp.gene.df)
      # Store the results in our list
      var_explained_summary_list[[result_index]] <- tibble(
        annotation = cur.annot,
        variant_id = cur.var,
        phenotype_id = cur.gene.name,
        rsquared = cur.rsq,
        pvalue = cur.pvalue,
        mean_expression = cur.mean.exp,
        median_expression = cur.median.exp,
        variance_expression = cur.var.exp,
        beta = cur.beta,
        beta.se = cur.beta.error,
        n_individuals = nrow(cur.gene.df),
        loci_ID = cur_eqtl_of_interest$loci_ID[1],
        Signal = cur_eqtl_of_interest$Signal[1],
        gwas_trait = cur_eqtl_of_interest$gwas_trait[1],
        Phenotype = cur_eqtl_of_interest$phenotype[1]
      )
      result_index <- result_index + 1
    }
  }
}

# Combine all results into one data frame
var.explained.df <- bind_rows(var_explained_summary_list)

if (SAVE.FILES == TRUE){
  write_tsv(var.explained.df,paste0(out.dir,"colocs_table-var_explained-",
                            gsub('[.]', 'pt', PP4.loose.threshold),".tsv"))
}
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

plot.var.explained.df <- var.explained.df %>% 
  annotate_sumstats()

# Make a dataframe for the top three cell types according to their rsquared at each locus
var.explained.top3.df <- var.explained.df %>% 
  annotate_sumstats() %>%
  group_by(variant_id) %>%
  slice_max(rsquared, n=3) %>% 
  # Make a category_new_count column based on the category_new with the max rsquared in the group
  mutate(category_new_count = paste0(category_new[which.max(rsquared)])) %>% 
  ungroup() %>% 
  group_by(category_new_count) %>%
  mutate(category_new_count = paste0(category_new_count, "\nn = ", n()/3)) %>% 
  ungroup() %>% 
  arrange(-rsquared) %>% 
  mutate(annotation = factor(annotation, levels=unique(annotation)))

# Get the IBD GWAS sumstats for aligning purposes
gwas.snps <- read_ibd_gwas('/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/gwas_sumstats_final/')


plot.gene <- 'PSEN2'
plot.gene <- 'ZMIZ1'
plot.gene.list <- c('MAML2', 'ZMIZ1', 'PSEN2', #Notch
                    'RASGRP1', 'LPIN3','MYC', 'FUBP1', 'FERMT1',# Wnt
                    'NDUFAF1', 'PRKCB', # Novel drug targets
                    'ITGA4', 'JAK2', 'ITGAL', 'TNFSF15', 'TNFRSF14') # Existing drug targets


# Base boxplots and variance explained
for (plot.gene in plot.gene.list){
  ##########
  # DATA WRANGLING
  ##########
  print(plot.gene)
  # For plot 1
  
  cur.coloc.df <- known.coloc.df %>% 
    filter(symbol == plot.gene, gwas_trait %in% ibd.traits) %>% 
    slice_max(PP.H4.abf)
  
  cur.plot.var.explained.df <- plot.var.explained.df %>% 
    filter(symbol == plot.gene, pvalue < 0.05, median_expression != 0) %>% 
    left_join(known.coloc.df %>% 
                filter(symbol == plot.gene,
                       gwas_trait == cur.coloc.df$gwas_trait[1],
                       PP.H4.abf > PP4.loose.threshold) %>%
                select(annotation) %>% 
                mutate(colocalising_annotation = TRUE),
                by = 'annotation')
  
  cur.var.explained.top.df <- cur.plot.var.explained.df %>% 
    slice_max(rsquared, n=3) %>%
    mutate(label_tissue = paste(label_new, tissue, sep = '\n'))
  
  # Get the GWAS SNP
  cur.gwas.snp <- gwas.snps %>% 
    filter(RSid == cur.var.explained.top.df$variant_id[1], 
           Disease == cur.coloc.df$gwas_trait[1])
  
  
  # For plot 2
  
  # Filter expression.paths for those that contain an annotation
  # from cur.var.explained.top.df$annotation
  cur.cell.types.paths <- int.expression.paths[str_detect(int.expression.paths, 
                                                          regex(paste0(cur.var.explained.top.df$annotation,
                                                                       collapse = '|')))]
  normalised <- TRUE
  
  cur.cell.types.dfs <- list()
  for (path in cur.cell.types.paths){
    # Load expression data 
    if (normalised){
      # Read dataframe and set names as a column
      cur.df <- read.table(path, header = TRUE, row.names = NULL) %>% 
        dplyr::rename(ensgene = row.names)
    } else {
      cur.df <- read.table(path,header = TRUE) %>% 
        dplyr::rename(ensgene = ENS)
    }
    cur.df <- cur.df %>% 
      filter(ensgene == cur.var.explained.top.df$phenotype_id[1]) %>%
      pivot_longer(-ensgene, names_to = 'annotation_individual') %>%
      separate(annotation_individual,
               into = c('annotation_type','annotation','individual'),
               sep = '[.]') %>% 
      mutate(part2 = str_split(individual, "_", simplify = TRUE)[,2],
             part3 = str_split(individual, "_", simplify = TRUE)[,3],
             genotype_id = if_else(part2 == part3, part2, paste(part2, part3, sep = "_"))) %>% 
      select(-part2, -part3, -individual)
    cur.cell.types.dfs[[path]] <- cur.df
  }
  cur.cell.types.df <- bind_rows(cur.cell.types.dfs)
  
  # Gene relevant SNP
  cur.snp <- vcf_lookup[[cur.var.explained.top.df$variant_id[1]]] %>%
    pivot_longer(-c(1:9), names_to = 'genotype_id') %>%
    mutate(genotype_id = str_remove(genotype_id, '^X')) %>%
    separate(value, into = c('GT','DS','HDS','GP'), sep = ':') %>%
    mutate(Genotype = str_count(GT, "1")) %>%
    select(ID, genotype_id, Genotype, GT) %>%
    drop_na()
  
  # Add genotype to the cur.cell.types.df
  plot.cur.cell.types.df <- left_join(cur.cell.types.df, cur.snp, by = 'genotype_id') %>% 
    dplyr::rename(phenotype_id = ensgene) %>% 
    mutate(variant_id = cur.var.explained.top.df$variant_id[1]) %>% 
    annotate_sumstats() %>% 
    mutate(label_tissue = factor(paste(label_new, tissue, sep = '\n'),
                                levels=rev(unique(cur.var.explained.top.df$label_tissue))))
  
  # Add doe adjusted to the gwas reference
  if (cur.gwas.snp$beta < 0){
    cur.plot.var.explained.df <- cur.plot.var.explained.df %>% 
      mutate(beta = beta * -1)
    
    plot.cur.cell.types.df <- plot.cur.cell.types.df %>%
      mutate(Genotype = abs(Genotype - 2))
  }
  
  ######
  # PLOT
  ######
  
  # Plot the rsqquared against the median expression
  p1 <- ggplot(cur.plot.var.explained.df,
               aes(x=rsquared, y=median_expression)) + 
    geom_pointrange(aes(ymin=median_expression-variance_expression,
                        ymax=median_expression+variance_expression,
                        fill=category_new,color=category_new,
                        shape=as.factor(sign(beta)))) +
    # stat_cor(label.x.npc = "center",label.y.npc = "top",colour='black') +
    geom_point(data=filter(cur.plot.var.explained.df, colocalising_annotation),
               aes(shape=as.factor(sign(beta))), colour='black',size=3.5,
               stroke = 1.5, show.legend = FALSE) +
    geom_text_repel(data=cur.var.explained.top.df,
                    aes(label=label_tissue,color=category_new),show.legend = FALSE,
                    size=2.5) +
    scale_x_continuous(labels = scales::percent, 
                       breaks = seq(0, 1, by = 0.05)) +
    scale_colour_manual(values = umap.category.palette,
                        name = 'Major\npopulation') +
    scale_fill_manual(values = umap.category.palette,
                      name = 'Major\npopulation') +
    scale_shape_manual(name = 'Effect of genetic variant\non gene expression',
                       values = c('1' = 24, 
                                  '-1' = 25),
                        labels = c("1"  = "Increase",
                                   "-1" = "Decrease")) +
    guides(shape = guide_legend(order = 1,nrow = 2),
           color = guide_legend(order = 2),
           fill = guide_legend(order = 2)) +
    labs(x='Variance in gene expression\nexplained by genotype', 
         y=paste0('Pseudobulked ', plot.gene,' expression')) +
    theme_classic() 
  
  # Make a boxplot faceted by cell type on the x-axis
  p2 <- ggplot(plot.cur.cell.types.df, aes(x=Genotype, y=value, 
                                               fill=category_new, group=Genotype)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(width=0.2, alpha=0.75,shape=16) +
    scale_fill_manual(values = umap.category.palette,
                      name = 'Major\npopulation',guide='none') +
    facet_wrap(~label_tissue, scales='free_x') +
    theme_classic() +
    scale_x_continuous(breaks = c(0,1,2), labels = c('0','1','2')) +
    labs(x=paste('Number of',cur.coloc.df$gwas_trait[1], 'risk increasing alleles'),
                 y=paste0('Rank-normalised\npseudobulked ',plot.gene,' expression'))
  
  
  # Combine plots and merge legends
  p.final <- ((p1 + p2) +
          plot_layout(guides = 'collect',widths = c(1,2)) & 
          plot_annotation(title = paste('Effect of',
                                        cur.var.explained.top.df$variant_id[1],
                                        'on', plot.gene, 'expression')) &
          theme(legend.position = 'bottom',legend.direction = 'horizontal'))
 
  
  # Save the plot
  if (SAVE.PLOTS == TRUE){
    ggsave(plot = p.final,paste0(out.dir,plot.gene,'_var_explained.pdf'),
           width = 11.6, height = 5,device = cairo_pdf)
  }
  
}


       