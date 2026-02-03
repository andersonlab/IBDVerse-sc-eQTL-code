library(tidyverse)
library(grid)

######################
# VARIABLES
######################

ibd.traits <- c('UC', 'CD', 'IBD')

level.mapping <- c('0'='All Cells', 
                   '1'='Major population',
                   '2'='Cell type')

tissue.mapping <- c("ct" = "Cross-site",
                    "ti" = "Terminal ileum" ,
                    "r" = "Rectum",
                    "blood" = "Blood")

interaction.mapping <- c('age_binned' = 'Age',
                         'disease_status' = 'Disease',
                         'ses_cd_binned' = 'Inflammation',
                         'ses_inflamed' = 'Inflammation2',
                         'sex' = 'Sex',
                         'smoking_ever' = 'Smoking')

PP4.loose.threshold <- 0.75
PP4.strict.threshold <- 0.9

######################
# DATA
######################


grch38 <- read_csv(paste0(repo.dir,'/data/grch38_genes.csv'))


annot.mapping <- read_csv(paste0(repo.dir,'/data/all_IBDverse_annotation_mastersheet.csv'), n_max = 103) %>% 
  dplyr::rename(label_machine = leiden,
                label_new = JAMBOREE_ANNOTATION,
                category_new = Category) %>% 
  dplyr::select(label_machine, label_new, category_new) %>%
  mutate(label_new = str_replace_all(label_new, '_', ' ')) %>% 
  mutate(Level = 2) %>% 
  bind_rows(., distinct(., category_new) 
            %>% mutate(label_new = category_new,
                       label_machine = category_new,
                       Level=1)) %>% 
  drop_na() %>% 
  add_row(label_machine='unannotated', label_new=unname(level.mapping['0']),
          category_new=unname(level.mapping['0']), Level=0) %>% 
  crossing(tissue = names(tissue.mapping)) %>% 
  mutate(label_machine = paste(label_machine, tissue, sep ='_'),
         tissue = tissue.mapping[tissue],
         annotation_type = level.mapping[as.character(Level)])

ibd.regions <- read_tsv(paste0(repo.dir,'/data/IBD_regions.txt')) %>% mutate(region_file='IBD')
cd.regions <- read_tsv(paste0(repo.dir,'/data/CD_regions.txt')) %>% mutate(region_file='CD')
uc.regions <- read_tsv(paste0(repo.dir,'/data/UC_regions.txt')) %>% mutate(region_file='UC')

known.ibd.regions <- bind_rows(ibd.regions, uc.regions, cd.regions) %>% 
  mutate(ibd_known_chr = as.numeric(chr),
         ibd_known_start = as.numeric(start),
         ibd_known_end = as.numeric(end),
         ibd_known_range = as.numeric(size)) %>% 
  dplyr::select(ibd_known_chr, ibd_known_start,ibd_known_end,
                ibd_known_range, loci_ID, Signal, Phenotype, region_file)

ibd.studies <- read_tsv(paste0(data.dir,'EFO_0003767_studies_export.tsv'))$accessionId

######################
# FUNCTIONS
######################

annotate_sumstats <- function(df){
  df <- df %>% 
  left_join(dplyr::select(grch38, ensgene, symbol, biotype),
            by = c("phenotype_id" = "ensgene"),multiple='first') %>% 
    mutate(symbol = case_when(symbol == '' ~ phenotype_id,
                              is.na(symbol) ~ phenotype_id,
                              TRUE ~ symbol)) %>% 
    mutate(variant_id = str_replace_all(variant_id, '_', ':')) %>% 
    separate(variant_id, into=c('chromosome', 'position', 'ref', 'alt'),
             sep=':',convert = TRUE,remove = FALSE) %>% 
    arrange(chromosome, position) %>% 
    left_join(annot.mapping, by=c('annotation' = 'label_machine'))
  return(df)
}

read_eqtls <- function(base_dir){
  sumstat.files <- list.files(base_dir, pattern="*Cis_eqtls_qval.tsv",
                              full.names=FALSE, recursive=TRUE)
  
  list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    name <- split.file.name[1]
    
    # Get info from name
    split.name <- str_split(name,'__')[[1]]
    annotation <- str_remove(split.name[2], '_all')
    method <- split.name[1]
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name))
    
    # Add new columns
    mutate(df,
           name = name,
           annotation = annotation,
           method = method,
           bonf_pvalue = pval_beta * n())
  })
  
  sumstat.df <- bind_rows(list.of.sumstat.dfs) %>% 
    left_join(dplyr::select(grch38, ensgene, symbol, biotype),
              by = c("phenotype_id" = "ensgene"),multiple='first') %>% 
    mutate(symbol = case_when(symbol == '' ~ phenotype_id,
                              TRUE ~ symbol)) %>% 
    left_join(annot.mapping, by=c('annotation' = 'label_machine'))
  
  geno_pheno <- list.files(paste0(base_dir,'../metadata/'),
                           '*genotype_phenotype_mapping.tsv',
                           full.names = TRUE) %>%
    map_dfr(read_tsv)
  
  # How many individuals per annotation
  n.individuals.df <- geno_pheno %>% 
    group_by(Sample_Category) %>% 
    tally(name = 'n_indiv') %>% 
    ungroup() %>% 
    separate(Sample_Category, sep = '-', into = c('annotation_type_raw',
                                                  'annotation',
                                                  'method'))
  
  # How many cells per individual-cell type pseudobulk
  # Filter out n < 5 because those are dropped before eQTL mapping is performed
  cell.metadata <- read_csv(paste0(base_dir,'../metadata/eqtl_processed.obs.csv'))
  cells.per.indiv.annot <- tibble()
  
  for (annot in unique(n.individuals.df$annotation_type_raw)){
    print(annot)
    tmp.df <- cell.metadata %>% 
      group_by(Genotyping_ID, get(annot)) %>% 
      mutate(n_cells=n()) %>% 
      filter(n_cells>=5) %>% 
      dplyr::slice(1) %>% 
      ungroup() %>% 
      mutate(annotation_type_raw = annot) %>% 
      dplyr::rename(annotation = `get(annot)`)
    cells.per.indiv.annot <- rbind(cells.per.indiv.annot, tmp.df)
  }
  
  # Get the average number of cells per annotation
  average.cells.per.annot <- cells.per.indiv.annot %>% 
    group_by(annotation) %>% 
    summarise(mean_n_cells = mean(n_cells)) %>% 
    ungroup()
  
  # Count number of genes expressed in each annot-indiv combination
  phenotype_file <- list.files(paste0(base_dir,'../metadata/'),
                               '*phenotype_file.tsv',
                               full.names = TRUE) %>%
    map_dfc(read_tsv)
  
  count.genes.df <- phenotype_file %>%
    summarize_all(~sum(. != 0))
  
  average.genes.per.annot <- count.genes.df %>% 
    dplyr::select(-starts_with('ENS')) %>% 
    pivot_longer(everything(), names_to = "annotation_sample", values_to = "n_genes") %>% 
    separate(annotation_sample,into = c('annotation', 'sample'), sep='-dMean_') %>% 
    group_by(annotation) %>% 
    summarise(mean_n_genes = mean(n_genes)) %>% 
    ungroup() %>% 
    separate(annotation, into = c('annotation_type_raw', 'annotation'), sep='-') %>% 
    dplyr::select(-annotation_type_raw)
  
  sumstat.df <- sumstat.df %>% 
    left_join(n.individuals.df, by = join_by(annotation)) %>%
    left_join(average.cells.per.annot, by = join_by(annotation)) %>%
    left_join(average.genes.per.annot, by = join_by(annotation)) %>%
    # left_join(disease.per.annot, by = join_by(Level, annotation))
  
  return(sumstat.df)
}

read_ieqtls <- function(base_dir){
  sumstat.files <- list.files(path=base_dir,
                              pattern="*cis_inter1.cis_qtl_top_assoc.txt.gz",
                              full.names=FALSE, recursive=TRUE)
  
  list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    name <- split.file.name[1]
    interaction <- split.file.name[4]
    print(paste(name, interaction))
    
    # Get info from name
    split.name <- str_split(name,'__')[[1]]
    annotation <- str_remove(split.name[2], '_all')
    method <- split.name[1]
  
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name))
    
    # Add new columns
    mutate(df,name = name,
           annotation = annotation,
           method = method,
           interaction = interaction)
  })
  
  sumstat.df <- bind_rows(list.of.sumstat.dfs) %>% 
    annotate_sumstats()  %>% 
    mutate(variant_id = factor(variant_id, levels=unique(variant_id)),
           interaction_new = interaction.mapping[interaction])
  
  return(sumstat.df)
}

# For when a file was saved for each GSEA result
read_gsea_many <- function(base_dir){
  sumstat.files <- list.files(path=base_dir,
                              pattern="*-gsea_results.tsv.gz",
                              full.names=FALSE)
  
  list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
    
    # Get info from name
    split.name <- str_split(file.name,'__')[[1]]
    annotation <- str_remove(split.name[2], '_all')
    method <- split.name[1]
    interaction <- split.name[3]
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name))
    
    # Add new columns
    mutate(df,name=file.name,
           annotation = annotation,
           method = method,
           interaction = interaction)
  })
  
  sumstat.df <- bind_rows(list.of.sumstat.dfs) %>%
    left_join(annot.mapping, by=c('annotation' = 'label_machine'))
  return(sumstat.df)
}

# For when the gsea results are in a single file
read_gsea_single <- function(base_dir){
  sumstat.files <- list.files(path=base_dir,
                              pattern="*-gsea_results.tsv.gz",
                              full.names=FALSE)
  
  list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name)) %>% 
      separate(coef_value, into=c('interaction', 'annotation_type', 'annotation', 'tissue'), sep=',') %>% 
      select(-tissue, -annotation_type) %>% 
      mutate(interaction_new = interaction.mapping[interaction]) %>% 
      left_join(annot.mapping, by=c('annotation' = 'label_machine'))
    
    # Add new columns
    mutate(df)
  })
  
  sumstat.df <- bind_rows(list.of.sumstat.dfs)
  return(sumstat.df)
}

read_conditional_eqtls <- function(base_dir){
  cond.ind.files <- list.files(base_dir, pattern="*Cis_eqtls_independent.tsv",
                               full.names=FALSE, recursive=TRUE)
  
  list.of.cond.ind.dfs <- lapply(cond.ind.files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    name <- split.file.name[1]
    
    # Get info from name
    split.name <- str_split(name,'__')[[1]]
    annotation <- str_remove(split.name[2], '_all')
    method <- split.name[1]
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name))
    
    # Add new columns
    mutate(df,name = name,
           annotation = annotation,
           method = method)
  })
  
  # Bind the independent signal DFs into a single DF and replace rank with the
  # order as the backwards step in the regression can result in gaps in the 
  # rank
  cond.ind.df <- bind_rows(list.of.cond.ind.dfs) %>% 
    annotate_sumstats() %>% 
    mutate(direction = sign(slope)) %>% 
    group_by(annotation_type, name, phenotype_id) %>%
    arrange(rank) %>% mutate(rank=row_number()) %>% ungroup()
  
  return(cond.ind.df)
}

read_trans <- function(base_dir){
  sumstat.files <- list.files(base_dir, pattern="*trans-by-cis_bonf_fdr.tsv",
                              full.names=FALSE, recursive=TRUE)
  
  list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    name <- split.file.name[1]
    
    # Get info from name
    split.name <- str_split(name,'-')[[1]]
    annotation.type <- split.name[1]
    annotation <- split.name[2]
    method <- split.name[3]
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name))
    
    # Add new columns
    mutate(df,
           name = name,
           annotation_type = annotation.type,
           annotation = annotation,
           method = method)
  })
  
  sumstat.df <- bind_rows(list.of.sumstat.dfs) %>%
    annotate_sumstats()
  
  return(sumstat.df)
}

read_pc_optim <- function(base_dir){
  pc.optim.files <- list.files(path=base_dir,
                              pattern="*optimise_nPCs-FDR0pt05.txt",
                              full.names=FALSE, recursive=TRUE)
  list.of.pc.optim.dfs <- lapply(pc.optim.files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    name <- split.file.name[1]
    
    # Get info from name
    split.name <- str_split(name,'__')[[1]]
    annotation <- str_remove(split.name[2], '_all')
    method <- split.name[1]
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name),col_types = 'ddl')
    
    # Add new columns
    mutate(df,
           name = name,
           annotation = annotation,
           method = method)
  })
  optim.pc.df <- bind_rows(list.of.pc.optim.dfs) %>% 
    left_join(annot.mapping, by=c('annotation' = 'label_machine')) %>%
    arrange(desc(num_eGenes)) %>%
    mutate(label_new = factor(label_new, levels=unique(label_new)))
  
  return(optim.pc.df)
  
}

read_doe <- function(base_dir){
  doe_files <- list.files(path=base_dir,
                              pattern="*coloc_doe.txt",
                              full.names=FALSE, recursive=TRUE)
  list.of.doe.dfs <- lapply(doe_files, function(file.name) {
    # Get info from filename
    split.file.name <- str_split(file.name,'/')[[1]]
    gwas <- split.file.name[1]
    
    # Read file
    df <- read_tsv(paste0(base_dir, file.name))
    
    # Add new columns
    mutate(df,
           gwas_trait = gwas)
  })
  # Combine all the files into a single dataframe
  doe.df <- bind_rows(list.of.doe.dfs) %>% 
    separate(condition_name, into=c('method',
                                      'annotation'), sep='__',remove = FALSE) %>%
    mutate(annotation = str_remove(annotation, '_all')) %>%
    left_join(annot.mapping, by=c('annotation' = 'label_machine'))
  
  return(doe.df)
}

map_var_to_region <- function(df){
  df %>% left_join(known.ibd.regions, by = c("chr" = "ibd_known_chr", 
                                        'region_file' = 'region_file'),
              relationship = "many-to-many") %>%
    dplyr::select(-region_file) %>% 
    filter(gwas_lead_pos >= ibd_known_start & gwas_lead_pos <= ibd_known_end)
}

get_colocs <- function(coloc.dir, known_ibd_only=FALSE){
  all.coloc.files <- list.files(coloc.dir, '*[.]gz')
  
  all.coloc.df.list <- lapply(all.coloc.files, function(file.name){
    df <- read_tsv(paste0(coloc.dir,'/', file.name)) %>% 
      separate(condition_name, into=c('method',
                                      'annotation'), sep='__',remove = FALSE) %>% 
      mutate(annotation = str_remove(annotation, '_all')) %>% 
      left_join(dplyr::select(grch38, ensgene, symbol, biotype),
                by = c("phenotype_id" = "ensgene"),multiple='first') %>% 
      mutate(symbol = case_when(symbol == '' ~ phenotype_id,
                                is.na(symbol) ~ phenotype_id,
                                TRUE ~ symbol)) %>% 
      filter(!(chr == '6' & gwas_lead_pos >= 28510120 & gwas_lead_pos <= 33480577)) # Filter HLA
  })
  
  
  all.coloc.df <- bind_rows(all.coloc.df.list) %>% 
    left_join(annot.mapping, by=c('annotation' = 'label_machine'))
  
  # Filters for known IBD regions only
  if (known_ibd_only == TRUE){
    
    known.coloc.df <- all.coloc.df %>%
      mutate(region_file = case_when(gwas_trait == 'CD' ~ 'CD',
                                      gwas_trait == 'UC' ~ 'UC',
                                      TRUE ~ 'IBD')) 
    known.coloc.df <- map_var_to_region(known.coloc.df)
    
  return(known.coloc.df)
  } 
  else {
    return(all.coloc.df)
    }
  
}

read_otar <- function(repo.dir){
  read_tsv(paste0(repo.dir,'data/all_otg_ibd_colocs.txt')) %>% 
    separate(Lead_variant,into = c('chr',
                                   'gwas_lead_pos',
                                   'ref',
                                   'alt'),sep = '_',remove = FALSE) %>% 
    mutate(region_file = case_when(grepl('Ulcerative',Trait_reported) ~ 'UC',
                                   grepl('Crohn', Trait_reported) ~ 'CD',
                                   TRUE ~ 'IBD'),
           chr=as.integer(chr),
           gwas_lead_pos=as.integer(gwas_lead_pos)) %>% 
    map_var_to_region()
}

read_ibd_gwas <- function(gwas.base.path){
  full.trait.list <- c('CD' = 'CrohnsDisease_DeLange_NatGen2017.txt.gz',
                       'UC' = 'UlcerativeColitis_DeLange_NatGen2017.txt.gz',
                       'IBD' = 'InflammatoryBowelDisease_DeLange_NatGen2017.txt.gz')
  
  cd.gwas.snps <- read_tsv(paste0(gwas.base.path,full.trait.list['CD'])) 
  uc.gwas.snps <- read_tsv(paste0(gwas.base.path,full.trait.list['UC'])) 
  ibd.gwas.snps <- read_tsv(paste0(gwas.base.path,full.trait.list['IBD']))
  
  # Concatenate the dataframes
  gwas.snps <- rbind(cd.gwas.snps, uc.gwas.snps, ibd.gwas.snps)
  return(gwas.snps)
}

get_ld_clumps_per_level <- function(ld.clump.path, cond.eqtl.df, sig.ieqtl.df=NULL){
  
  if (is.null(sig.ieqtl.df)) {
    all.eqtls <- cond.eqtl.df
  } else {
    all.eqtls <- bind_rows(select(cond.eqtl.df, phenotype_id, variant_id, slope, label_new, Level),
                           select(rename(sig.ieqtl.df, slope = b_gi),
                                  phenotype_id, variant_id, slope, label_new, Level, interaction))
  }
  print(nrow(all.eqtls))
  ld_clumps <- read_table(ld.clump.path) %>% 
    inner_join(all.eqtls, by=c('variant_id', 'phenotype_id')) 
  
  print(nrow(ld_clumps))
  
  ld_indices_min_level <- ld_clumps %>% 
    group_by(qtl_clump_index, phenotype_id) %>%
    distinct(label_new,.keep_all = TRUE) %>%  # Collapses same eQTL, same annot in different tissues
    mutate(n_level0 = sum(Level == 0),
           n_level1 = sum(Level == 1),
           n_level2 = sum(Level == 2)) %>%
    slice_min(Level, with_ties = FALSE) %>%
    ungroup() %>% 
    separate(qtl_clump_index, into = c('chr', 'pos', 'ref', 'alt'), sep = ':',remove = FALSE) %>% 
    mutate(pos = as.integer(pos))
  
  print(nrow(ld_indices_min_level))
  return(ld_indices_min_level)
}

read_monogenic <- function(repo.dir){
  monogenic.ibd.genes <- read_csv(paste0(repo.dir,'data/list_monogenic_ibd_list.csv')) %>% 
    select(Gene, list, `Effect  (Uhlig JPGN 2020)`, `Inheritance AR/AD/XL  (Uhlig JPGN 2020)`) %>% 
    mutate(region_file = 'IBD') %>% 
    left_join(filter(grch38, chr %in% c(1:22, 'X')), by=c('Gene'='symbol')) %>% 
    filter(chr != 'X') %>% # Remove X as we don't use it
    mutate(gwas_lead_pos = start + (end - start)/2,
           chr = as.numeric(chr)) %>% # No actual gwas lead pos so using midpoint as proxy
    map_var_to_region()
}

# Helper function to get the label from whichever ontology is used
get_term_label_for_id <- function(ont_id) {
  if (stringr::str_detect(ont_id, "^MONDO")) {
    return(termLabel(Term(mondo, ont_id)))
  } else if (stringr::str_detect(ont_id, "^EFO")) {
    return(termLabel(Term(efo, ont_id)))
  } else if (stringr::str_detect(ont_id, "^Orphanet")) {
    return(termLabel(orphanet.terms@x[[ont_id]]))
  } else if (stringr::str_detect(ont_id, "^HP")) {
    return(termLabel(Term(hp, ont_id)))
  } else {
    return("Unknown Ontology")
  }
}

# Classify disease based on its ontology ID
classify_disease <- function(ontology_id) {
  # Ontology ancestral terms 
  # Aim of these is to have terms high up in the ontology that
  # the diseases can be broadly binned in to
  # Order of diseases matters here so an immune disease that is also nutritional
  # would be labelled immune, in current ordering
  ONTOLOGY_CATEGORIES <- list(
    # IBD is the only "specific disease"
    IBD = c('EFO:0003767',      # 'inflammatory bowel disease' in EFO
            'MONDO:0005265'),      # 'inflammatory bowel disease' in MONDO
    INFECTION = c("EFO:0005741",     # 'infectious disease' in EFO
                  "MONDO:0005550"),     # 'infectious disease' in MONDO
    NEOPLASTIC = c("MONDO:0004992",   # 'neoplastic disease' in MONDO
                   "MONDO:0045024",   # 'malignant neoplasm' in MONDO
                   "EFO:0000311"),     # 'cancer' in EFO
    IMMUNE_SYSTEM = c("EFO:0000540",      # 'immune system disease' in EFO
                      'MONDO:0005046'),     # 'immune system disorder' in MONDO
    METABOLIC = c("MONDO:0005066",   # 'metabolic disease' in MONDO
                  "EFO:0000589"),     # 'metabolic disease' in EFO
    NUTRITIONAL = c("EFO:0001069",     # 'nutritional disease' in EFO
                    "MONDO:0005137"),     # 'nutritional disorder' in MONDO
    CARDIOVASCULAR = c("EFO:0000319",     # 'cardiovascular disease' in EFO
                       "MONDO:0004995"),     # 'cardiovascular disorder' in MONDO
    ENDOCRINE = c("MONDO:0005151",   # 'endocrine system disease' in MONDO
                  "EFO:0001379"),     # 'endocrine system disease' in EFO
    RESPIRATORY = c("MONDO:0005087",   # 'respiratory system disorder' in MONDO
                    "EFO:0000684"),     # 'respiratory system disease' in EFO
    MUSCULOSKELETAL = c("MONDO:0002081",   # 'musculoskeletal system disease' in MONDO
                        "EFO:0009676"),     # 'musculoskeletal system disease' in EFO
    NEUROLOGICAL = c("MONDO:0005071",   # 'nervous system disorder' in MONDO
                     "EFO:0000618",      # 'nervous system disease' in EFO
                     "ORDO:98006", # 'rare neurologic disease' in ORDO
                     "Orphanet:98006"),     # 'rare neurologic disease' in ORDO
    DIGESTIVE = c("MONDO:0004335",   # 'digestive system disorder' in MONDO
                  "EFO:0000405"),     # 'digestive system disease' in EFO
    RENAL = c("EFO:0003086",    # 'kidney disease' in EFO
              "MONDO:0005240"),     # 'kidney disorder' in MONDO
    HEPATIC = c("MONDO:0005154",   # 'liver disorder' in MONDO
                'HP:0410042',
                "EFO:0001421"),     # 'liver disease' in EFO
    DERMATOLOGICAL = c("MONDO:0005093",   # 'skin disorder' in MONDO
                       "EFO:0000701"),     # 'skin disease' in EFO
    HEMATOLOGICAL = c("MONDO:0005570",   # 'hematologic disease' in MONDO
                      "EFO:0005803"),     # 'hematological disease' in EFO
    MENTAL_BEHAVIORAL = c("MONDO:0005084",   # 'mental disorder' in MONDO
                          "EFO:0000677",     # 'mental or behavioural disorder' in EFO
                          "HP:0011446"),     # 'Abnormality of mental function' in MONDO
    SUBSTANCE_USE = c("MONDO:0002491", # Substance abuse
                      "MONDO:0002494"), # Substance related disorder
    EYE = c('EFO:0003966', #'eye disease' in EFO
            'MONDO:0005328'), # 'eye disorder' in MONDO
    POISONING = c("MONDO:0029000",   # 'poisoning' in MONDO
                  "EFO:0008546"),     # 'poisoning' in EFO
    INJURY = c("MONDO:0021178",   # 'injury' in MONDO
               "EFO:0000546")     # 'injury' in EFO
  )
  print(ontology_id)
  # 1) Basic checks for empty / unknown IDs
  if (is.na(ontology_id) || !nzchar(ontology_id)) {
    return(c(ontology_id, "Other", "No ID"))
  }
  
  # 2) Decide which ontology to query
  if (str_detect(ontology_id, "^MONDO")) {
    term  <- Term(mondo, ontology_id)
    label <- termLabel(term)
    ancs  <- tryCatch(ancestors(term), error = function(e) NULL)
  } else if (str_detect(ontology_id, "^EFO")) {
    term  <- Term(efo, ontology_id)
    label <- termLabel(term)
    ancs  <- tryCatch(ancestors(term), error = function(e) NULL)
  } else if (str_detect(ontology_id, "^Orphanet")) {
    ontology_id <- str_replace_all(ontology_id, "_", ":")
    term <- orphanet.terms@x[[ontology_id]]
    if (is.null(term)){
    ontology_id <- str_replace_all(ontology_id, "Orphanet", "ORDO")
    term <- orphanet.terms@x[[ontology_id]]}
    label <- termLabel(term)
    ancs  <- tryCatch(ancestors(term), error = function(e) NULL)
  } else if (str_detect(ontology_id, "^HP")) {
    term  <- Term(hp, ontology_id)
    label <- termLabel(term)
    ancs  <- tryCatch(ancestors(term), error = function(e) NULL)
  } else {
    # Not MONDO/EFO
    return(c(ontology_id, "Other", "Not MONDO/EFO/HP/Orphanet"))
  }
  
  # If ancestors can't be fetched, fallback to Other
  if (is.null(ancs)) {
    return(c(ontology_id, "Other", label))
  }
  
  # 3) Combine the disease's own ID + all ancestor IDs
  ancs_ids <- map_chr(ancs, "obo_id", .default = NA_character_)
  all_ids  <- unique(c(ancs_ids, str_replace_all(ontology_id, '_', ':')))
  
  # 4) Loop through each category in CATEGORIES
  for (cat_name in names(ONTOLOGY_CATEGORIES)) {
    top_ids  <- ONTOLOGY_CATEGORIES[[cat_name]]
    
    # Check if any of that category's top-level IDs is in the disease's lineage
    if (any(top_ids %in% all_ids)) {
      # The user wants the second item in the return vector to be the label from 
      # the first ID in `top_ids`. We'll fetch that from the ontology:
      cat_label <- get_term_label_for_id(top_ids[1])
      # Return (ontology_id, categoryLabelFromOntology, diseaseLabel)
      return(c(ontology_id, cat_label, label))
    }
  }
  
  # 5) If we don't match any category => default to "Other"
  return(c(ontology_id, "Other", label))
}

# Function to plot coloc at a locus for all colocs in a given coloc dataframe
plot_coloc_locus <- function(trait, coloc_df,gwas_sumstats_path,eqtl_sumstats_base_path,
                             interaction.or.base,doe_or_plot,pph4_thresh,out_dir,
                             ld_indices_min, plot_gene_list=NULL,plot.window.size=5e5,
                             all_variants=FALSE){

  library(ggrepel)
  
  # Filter coloc df for a trait
  coloc_df <- coloc_df %>% filter(gwas_trait == trait, PP.H4.abf >= pph4_thresh)

  # Find highest PP4 per annot type
  gwas.df <- read_tsv(gwas_sumstats_path)
  
  if (nrow(coloc_df) == 0){
    print('No colocs with this GWAS')
  } else{ 
    for (i in 1:nrow(coloc_df)){
      print(i)
      # Subset for the row
      row <- coloc_df[i,]
      if (interaction.or.base == 'base'){
        if (row$interaction != 'base'){
          print('Not a base coloc')
          next
        }
      } else if (interaction.or.base == 'interaction'){
        if (row$interaction == 'base'){
          print('Not an interaction coloc')
          next
        }
      }
      
      # Extract relevant info from the row
      var.id <- row$gwas_lead
      condition_name <- row$condition_name
      ens.gene <- row$phenotype_id
      gene.id <- row$symbol
      PP4 <- round(row$PP.H4.abf,2)
      annotation <- row$annotation
      tissue <- row$tissue
      label_new <- row$label_new
      print(condition_name)
      print(gene.id)
      if (!is.null(plot_gene_list)){
        if (gene.id %in% plot_gene_list){
          print('Gene in list')
        } else {
          next
        }
      }
      
      # Variable to minimise dataframe loading
      prev.condition_name <- ''
      
      
      PP4.thresh.text <- gsub('[.]', 'pt', paste0('PPH4_',pph4_thresh))
      file.out.path <- paste(trait,
                             gene.id,
                            tissue,
                            var.id, 
                            condition_name,
                            gsub('[.]', 'pt', paste0('PPH4_',PP4)),
                       sep = '-', collapse = '')
      if (doe_or_plot == 'plot'){
        loci.out.dir <- paste0(out_dir,"coloc_loci/",
                               PP4.thresh.text,interaction.or.base,"-all_vars"[all_variants],'/')
        file.out.path <- paste0(loci.out.dir, '/',file.out.path,".pdf")
      } else {
        loci.out.dir <- paste0(out_dir,"coloc_doe/",PP4.thresh.text,interaction.or.base,"-all_vars"[all_variants],'/')
        file.out.path <- paste0(loci.out.dir, '/',file.out.path,".txt")
      }
      
      print(file.out.path)
      # Check if file already exists
      if (file.exists(file.out.path)) {
        print('Coloc file has already been made')
        # If it does, skip to the next iteration
        next
      }
      
      var.id.chromosome <- gsub('chr', '', str_split(var.id, '_')[[1]][1])
      var.id.pos <- as.numeric(str_split(var.id, '_')[[1]][2])
      
      if (interaction.or.base == 'base'){
        eqtl_sumstats_path <- paste0(eqtl_sumstats_base_path,
                                     condition_name,
                                     '/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.',var.id.chromosome,'.tsv')
      } else if (interaction.or.base == 'interaction'){
        interaction <- str_split(condition_name,pattern = '__')[[1]][3]
        condition_name <- sub("__[^_]+(_[^_]+)*$", "", condition_name)
        eqtl_sumstats_path <- paste0(eqtl_sumstats_base_path,
                                     condition_name,
                                     '/OPTIM_pcs/interaction_output/',
                                     interaction
                                     ,'/cis_inter1.cis_qtl_pairs.',var.id.chromosome,'.tsv')
      }
      
      eqtl.df <- read_tsv(eqtl_sumstats_path,show_col_types = FALSE)
      
      chr.eqtl.df <- eqtl.df %>%  
        separate(variant_id, c('Chr', 'Pos', 'Ref', 'Alt'),':') %>% 
        mutate(Chr = as.numeric(gsub('chr', '', Chr)),
               Pos = as.numeric(Pos)) %>% 
        filter(phenotype_id == ens.gene)
      
      chr.gwas.df <- gwas.df %>% 
        filter(Chr == var.id.chromosome)
      
      if (!all_variants){
        # Do a left join and drop_na so variants not in both dataframes are dropped
        joint.df <- left_join(chr.eqtl.df, gwas.df, by=c('Pos', 'Chr')) %>% 
          mutate(beta = case_when(Eff_allele == Ref ~ -beta,
                                  Eff_allele == Alt ~ beta,
                                  TRUE ~ NA)) %>% 
          drop_na(beta)
      } else {
        # Keep variants that are in either dataframe
        joint.df <- full_join(chr.eqtl.df, gwas.df, by=c('Pos', 'Chr')) %>% 
          mutate(beta = case_when(Eff_allele == Ref ~ -beta,
                                  Eff_allele == Alt ~ beta,
                                  TRUE ~ NA))
      }
      
      min_bp <- var.id.pos - plot.window.size
      max_bp <- var.id.pos + plot.window.size
      plot.joint.df <- filter(joint.df, (Pos > min_bp) &  (Pos < max_bp))
      
      # Make a dataframe for labelling the coloc genes
      label_colocs <- coloc_df %>% 
        group_by(phenotype_id) %>% 
        slice_max(PP.H4.abf, with_ties = FALSE) %>% 
        select(phenotype_id, category_new)
      
      # Plot all the genes in the locus
      gene.track.df <- filter(grch38, (end > min_bp) &  (start < max_bp),
                              chr == var.id.chromosome) %>% 
        arrange(start) %>%               # sort by start
        { # We keep a local 'track_ends' vector that holds, for each track,
          # the 'end' coordinate of the last gene placed on that track.
          track_ends <- numeric(0)       # initially empty (no tracks)
          # We use map2_int (from purrr) to iterate over (start,end) pairs
          mutate(., track = map2_int(start, end, ~ {
                   # Find first track whose 'end' is < this gene's 'start'
                   idx <- which(.x > track_ends)[1]
                   if (is.na(idx)) {
                     # no existing track is free, so we add a new one
                     track_ends[length(track_ends) + 1] <<- .y
                     length(track_ends)  # this new track index
                   } else {
                     # reuse track idx
                     track_ends[idx] <<- .y
                     idx }})) } 
      gene.track.df <- gene.track.df %>% 
        left_join(label_colocs,by = c('ensgene'='phenotype_id')) %>% 
        mutate(is_coloc = case_when(is.na(category_new) ~ FALSE,
                                    TRUE ~ TRUE),
               # Dependent on strand, swap the start and end
               start = ifelse(strand == '1', start, end),
               end = ifelse(strand == '1', end, start))
      
      plot.joint.df$loggwaspvalue <- -log10(plot.joint.df$pval)
      if (interaction.or.base == 'interaction'){
        plot.joint.df$logeqtlpvalue <- -log10(plot.joint.df$pval_gi)
        plot.joint.df$slope <- plot.joint.df$b_gi
      } else{
        plot.joint.df$logeqtlpvalue <- -log10(plot.joint.df$pval_nominal)
      }
      pval_scaler = max(plot.joint.df$loggwaspvalue, na.rm = TRUE) / max(plot.joint.df$logeqtlpvalue,na.rm = TRUE)
      
      # Get GWAS lead variant
      lead.var.df <- plot.joint.df  %>% 
        drop_na(loggwaspvalue, logeqtlpvalue) %>%
        slice_max(logeqtlpvalue, n=1)
      
      if (doe_or_plot == 'plot'){
        # Colours
        gwas_colour <- 'darkgrey'
        eqtl_colour <- umap.category.palette[[row$category_new]]
        # eqtl_colour <- 'black'
        
        
        ld_indices_min_level_gene <- ld_indices_min %>% 
          filter(phenotype_id == ens.gene) %>% 
          mutate(pos = as.numeric(pos),
                 lvl2_clamped = factor(pmin(n_level2, 5), levels = 1:5),
                 lvl1_clamped = factor(pmin(n_level1, 5), levels = 1:5),
                 lvl0_clamped = factor(pmin(n_level0, 5), levels = 1:5))
        
        dummy_ld_indices_min_level_gene <- data.frame(
               pos          = NA_real_,
               lvl2_clamped = factor(as.character(1:5), levels = as.character(1:5)),
               y = unname(level.mapping['2']))
        
        p1 <- ggplot(ld_indices_min_level_gene, aes(x=pos)) +
          geom_point(data = dummy_ld_indices_min_level_gene,
                     aes(x = pos, y = y, colour = lvl2_clamped),
            shape = 124, size = 4, show.legend= TRUE) + # Dummy so legend works
          geom_point(data=filter(ld_indices_min_level_gene, n_level2 > 0),
                     aes(y=unname(level.mapping['2']), colour=lvl2_clamped), size=4, shape=124) +
          geom_point(data=filter(ld_indices_min_level_gene, n_level1 > 0),
                     aes(y=unname(level.mapping['1']), colour=lvl1_clamped), size=4, shape=124) +
          geom_point(data=filter(ld_indices_min_level_gene, n_level0 > 0),
                     aes(y=unname(level.mapping['0']), colour=lvl0_clamped), size=4, shape=124) +
          scale_colour_viridis_d(
            name   = paste0('Number of\nannotations\nwith an independent\n', gene.id, ' eQTL'),
            limits = as.character(1:5), breaks = as.character(1:5), labels = c("1","2","3","4","5+"),
            drop   = FALSE, direction = -1, end = 0.8) +
          scale_shape_manual(
            name   = paste0('Number of\nannotations\nwith an independent\n', gene.id, ' eQTL'),
            values = rep(124, 5), limits = as.character(1:5), labels = c("1","2","3","4","5+"),
            drop   = FALSE) +
          scale_y_discrete(limits=c(level.mapping['2'], level.mapping['1'], level.mapping['0'])) +
          scale_x_continuous(limits = c(min_bp, max_bp)) +
          theme_classic() +
          guides(colour=guide_legend(ncol=3, 
                                     override.aes = list(shape = 124, size = 4))) +
          theme(axis.title = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.box = 'horizontal',
                legend.key.size = unit(0.25, 'cm'),
                legend.title = element_text(size=6),
                legend.text = element_text(size=5),
                plot.margin = unit(c(0,0,0,0), "cm"))  +
          labs(title = paste0(gene.id,", PPH4 = ", PP4))
       
        # Titles for p2
        p2_eqtl_y_title <- paste0("-log[10](P[eQTL - ", shQuote(label_new), "])")
        p2_gwas_y_title <- paste0("-log[10](P[GWAS - ", shQuote(trait), "])")
        
        p2 <- ggplot(plot.joint.df) +
          geom_point(aes(x=Pos, y=loggwaspvalue, color="GWAS"),size=3)  +
          geom_point(aes(x=Pos, y=logeqtlpvalue*pval_scaler, colour="eQTL"),size=3,alpha=0.5)  +
          scale_color_manual(values=c("GWAS"=gwas_colour, "eQTL"=eqtl_colour), name=paste0("PPH4 = ", PP4)) +
          scale_y_continuous(sec.axis = sec_axis(~./pval_scaler, name=parse(text=p2_eqtl_y_title)),
                             name = parse(text=p2_gwas_y_title)) +
          scale_x_continuous(labels = function(x) paste0(round(x / 1e6, 1), " Mb"), 
                             name = paste('Chromosome',var.id.chromosome),
                             limits = c(min_bp, max_bp)) +
          theme_classic() +
          theme(legend.position = 'none',
                axis.title.y.left = element_text(colour = gwas_colour,face = 'bold',size=10),
                axis.title.y.right = element_text(color = eqtl_colour,face = 'bold',size=10),
                axis.text.y.left = element_text(colour = gwas_colour,size=8),
                axis.text.y.right = element_text(color = eqtl_colour,size=8),
                axis.title.x = element_blank(), 
                axis.text.x = element_blank(),
                axis.ticks.x=element_blank())
        
        p3 <- ggplot(gene.track.df, aes(y= track,x= start,
                                        xend = end,yend = track,
                                        colour=is_coloc)) +
          geom_segment(arrow = arrow(length=unit(0.1, "cm")), size=0.5) +
          geom_text_repel(aes(x = (start + end)/2, label = symbol,y=track),
                          size = 2,fontface = "bold") +
          scale_colour_manual(values=c('TRUE' = 'black', 
                                       'FALSE' = 'grey'),
                              guide='none') + 
          scale_x_continuous( labels = function(x) paste0(round(x / 1e6, 1), " Mb"),
                              name = paste('Genomic position on chromosome',var.id.chromosome),
                              limits = c(min_bp, max_bp)) +
          scale_y_continuous(limits = c(0, 1+max(gene.track.df$track)),) +
          theme_classic() + 
          theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                legend.position = 'None')

        
        # p1/p2/p3 + plot_layout(heights = c(4,15,5))
        # # needs patchwork 1.3: library(patchwork, lib.loc = '/software/team152/oa3/r_libs/cran/4.3.1/')
        p.final <- (p1/free(p2, 'label')/p3 + plot_layout(heights = c(4,15,3)))
        
        ggsave(plot = p.final,file.out.path,width = 8, height = 5,device = cairo_pdf)
      } else if (doe_or_plot == 'doe') {
        # Calculate the doe based on lead variant
        lead.var.beta.ratio <- sign(lead.var.df$beta /lead.var.df$slope)
        # Calculate effect based on average across all variants in locus where
        # beta is present for eqtl AND gwas
        doe.df <- drop_na(plot.joint.df, loggwaspvalue, logeqtlpvalue) %>% 
          filter(Pos > var.id.pos - 1e5, Pos < var.id.pos + 1e5,
                 loggwaspvalue > 2, logeqtlpvalue > 2)
        all.var.doe <- sign(doe.df$beta / doe.df$slope)
        average.beta.ratio.sign <- mean(all.var.doe)
        # Write out the file
        tibble(condition_name = condition_name,
               gwas=trait, gene = gene.id,
               lead_var_doe = lead.var.beta.ratio, 
               locus_doe = average.beta.ratio.sign,
               n_variants=nrow(doe.df)) %>% 
          write_tsv(file.out.path)
      }
      
      prev.condition_name <- condition_name
    }
  }
}


# Function to wrap labels for plotting
wrap_label <- function(x) str_wrap(str_replace_all(x, "_", " "), width = 25)

# Function to split up long strings
insert_newline_at_midpoint <- function(text, min_length = 50) {
  # Only process if the string length is at least min_length.
  if(nchar(text) < min_length) return(text)
  
  # Find the positions of all spaces in the string.
  space_positions <- gregexpr(" ", text)[[1]]
  
  # If no spaces found, simply return the original text.
  if(space_positions[1] == -1) return(text)
  
  # Calculate the midpoint of the string.
  mid <- nchar(text) / 2
  
  # Identify which space is closest to the midpoint.
  idx <- which.min(abs(space_positions - mid))
  pos <- space_positions[idx]
  
  # Insert a newline in place of that space.
  new_text <- paste0(substr(text, 1, pos - 1), "\n", substr(text, pos + 1, nchar(text)))
  new_text
}

###############
# PALETTES
###############

# Set up the colour palette for unnannotated/bulk-like, category/major population
# and cell type
annot.class.palette <- setNames(
  c('#edf8b1', '#7fcdbb', '#2c7fb8'),
  level.mapping
)

# Set up the point size for unannotated/bulk-like, category/major population
# and cell type
annot.class.sizes <- setNames(
  c(12, 7, 3),
  level.mapping
)

umap.palette.df <- read_csv(paste0(repo.dir,'/data/palette.csv'))
umap.category.palette <- c(deframe(dplyr::select(umap.palette.df, category, category_color)), 
                           annot.class.palette)
umap.celltype.palette <- c(deframe(dplyr::select(umap.palette.df, manual_annotation, manual_annotation_color)), 
                           annot.class.palette)

inflammation.palette <- c(uninflamed='#ffffbf',mild='#fdae61',
                          moderate='#f46d43',severe='#d73027')

disease.palette <- c('Healthy' = '#353D6D', 'CD' = '#F7B817')


otar.palette <- c('FALSE'='#3489ca', 'TRUE'='#acd1e8')

tissue.palette <- list(
  "Blood"                          = "#B15663",    # Solid
  "Rectum"                         = "#E07C2E",    # Solid
  "Terminal ileum"                 = "#30844D",    # Solid
  "Cross-site"                      = "grey"       # Solid
)

tissue.pattern.palette <- list(
  # Two-way intersections: linear blend of the two tissue colors
  "Blood, Rectum"                  = linearGradient(c(tissue.palette[['Blood']], 
                                                      tissue.palette[['Rectum']]),
                                                    stops=c(0.05,0.08)),
  "Rectum, Terminal ileum"         = linearGradient(c(tissue.palette[['Rectum']], 
                                                      tissue.palette[['Terminal ileum']]),
                                                    stops=c(0.5,0.8)),
  "Blood, Terminal ileum"          = linearGradient(c(tissue.palette[['Blood']], 
                                                      tissue.palette[['Terminal ileum']]),
                                                    stops=c(0.15,0.25)),
  
  # Three-way intersection: radial blend of all three tissue colors
  "Blood, Rectum, Terminal ileum"  = linearGradient(c(tissue.palette[['Blood']],
                                                      tissue.palette[['Rectum']],
                                                      tissue.palette[['Terminal ileum']]))
)

tissue.palette <- c(tissue.palette,tissue.pattern.palette)

