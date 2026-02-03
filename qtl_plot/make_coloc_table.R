library(tidyverse)
library(patchwork)
library(xtable)

##############
# Variables
###############

repo.dir <- '~/eqtl/code/IBDVerse-sc-eQTL-code/'

coloc.base.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/2025_03_IBDverse_coloc_all_gwas/collapsed/'
coloc.interaction.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed'

out.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/'



data.dir <- paste0(repo.dir,'/data/')
# SAVE
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE


##############
# Functions
##############

source(paste0(repo.dir,'/qtl_plot/helper_functions.R'))

##############
# Get the data
###############

known.base.coloc.df <- get_colocs(coloc.base.dir, known_ibd_only = TRUE)
known.interaction.coloc.df <- get_colocs(coloc.interaction.dir, known_ibd_only = TRUE) %>% 
  mutate(interaction = str_split_fixed(condition_name, "__", 3)[,3])

known.coloc.df <- bind_rows(mutate(known.base.coloc.df, interaction = 'base'),
                                        known.interaction.coloc.df)

otar.ibd.colocs <- read_otar(repo.dir)

otar.hits <- unique(filter(otar.ibd.colocs,H4 > PP4.loose.threshold)$Gene_symbol)

monogenic.ibd.genes <- read_monogenic(repo.dir)

####################
# All coloc/novel hits
####################

# Make a table with a row per locus and each row contains the following columns
# The GWAS trait, the locus ID, the lead SNP, the max PPH4 for the locus, 
# all colocalising genes (concatenated with a comma separating),
# the number of genes in the locus, whether the locus has an OTAR genetics hit,
# Whether any of the colocalising genes are in the monogenic IBD list of genes
# all of the colocalising genes at the same locus when using the new gwas (ibd, uc, cd)

make_coloc_table <- function(sig.ibd.coloc.hits){
  table.coloc.hits.df <- sig.ibd.coloc.hits %>%
    select(gwas_trait, gwas_lead, loci_ID, PP.H4.abf, symbol, interaction) %>%
    group_by(loci_ID, symbol) %>%
    mutate(max_PPH4 = max(PP.H4.abf)) %>% 
    ungroup() %>% 
    arrange(desc(PP.H4.abf)) %>%
    group_by(loci_ID) %>%
    mutate(monogenic_ibd_genes = ifelse(symbol %in% monogenic.ibd.genes$Gene, symbol, NA),
           OTG_ibd_coloc_genes = ifelse(symbol %in% otar.ibd.colocs$Molecular_trait, symbol, NA)) %>% 
    select(-PP.H4.abf) %>%
    summarise(across(
      .fns  = ~ paste(unique(na.omit(.)), collapse = ", ")),
      .groups = "drop") %>% 
    ungroup() %>% 
    # Round all occurrences of floating point numbers in the string column max_PPH4
    mutate(
      max_PPH4 = str_split(max_PPH4, ",\\s*") %>% 
        map_chr(function(x) {
          nums    <- as.numeric(x)
          rounded <- signif(nums, 3)
          # format to exactly 3 decimal places
          formatted <- formatC(rounded, format = "f", digits = 3)
          paste(formatted, collapse = ", ")
        })
    )
}

delangemoutsilee.coloc.hits <- make_coloc_table(filter(known.coloc.df,
                                              gwas_trait %in% ibd.traits,
                                              PP.H4.abf > PP4.loose.threshold))

fachal.coloc.hits <- make_coloc_table(filter(known.coloc.df,
                                        gwas_trait %in% c('uc', 'cd', 'ibd'),
                                        PP.H4.abf > PP4.loose.threshold))

table.coloc.hits.df <- delangemoutsilee.coloc.hits %>% left_join(select(fachal.coloc.hits, loci_ID, symbol),
                                 by = c('loci_ID'), suffix = c('.delange', '.fachal'))
  


if (SAVE.FILES == TRUE){
  write_tsv(delangemoutsilee.coloc.hits,paste0(out.dir,"colocs_table-PPH4gt",
                                       gsub('[.]', 'pt', PP4.loose.threshold),".tsv"))
}


# Filter for the novel colocs (not in OTG or monogenic list) with high confidence
# and where only a single effector gene is nominated at the locus and then 
# make a LaTeX table of the results
delangemoutsilee.coloc.hits %>% 
  filter(monogenic_ibd_genes == "",
         OTG_ibd_coloc_genes == "",
         max_PPH4 > PP4.strict.threshold,
         !str_detect(symbol, ',|ENSG')) %>%
  mutate(max_PPH4 = as.numeric(max_PPH4),
         loci_ID = as.character(loci_ID)) %>% 
  select(-OTG_ibd_coloc_genes,-monogenic_ibd_genes) %>%
  xtable()
  

####################
# Variance explained per coloc
####################

# Load in the variance explained dataframe
variance.explained.df <- read_tsv(paste0(out.dir, 'coloc_loci/colocs_table-var_explained-0pt75.tsv')) %>% 
  annotate_sumstats()

# Variance explained interpretable table
table.variance.explained <- variance.explained.df %>% 
  filter(rsquared > 0.05, pvalue < 0.05) %>% 
  mutate(variance_explained = round(rsquared, 3),
         median_expression = round(median_expression, 3)) %>%
  select(tissue, label_new, category_new, loci_ID, symbol, variance_explained, median_expression) %>%
  arrange(desc(variance_explained)) %>%
  group_by(loci_ID, symbol) %>% 
  slice_max(variance_explained, n=3) %>%
  summarise(across(
    .fns  = ~ paste(na.omit(.), collapse = ", ")),
    .groups = "drop") %>% 
  ungroup()


if (SAVE.FILES == TRUE){
  write_tsv(table.variance.explained,paste0(out.dir,"variance_explained_table-PPH4gt",
                                       gsub('[.]', 'pt', PP4.loose.threshold),".tsv"))
}
