library(tidyverse)
library(otargen)

##############
# Variables
###############

# server <- 'farm'
server <- 'openstack'

if (server == 'farm'){
  repo.dir <- '~/eqtl/code/IBDVerse-sc-eQTL-code/'
  gwas.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/snakemake_colocalisation/gwas_sumstats/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  gwas.dir <- ''
}


data.dir <- paste0(repo.dir,'/data/')
# SAVE
SAVE.PLOTS <- TRUE
SAVE.FILES <- TRUE

MAKE.REGIONS <- FALSE


##############
# Functions
##############

source(paste0(repo.dir,'/qtl_plot/helper_functions.R'))

create_non_overlapping_regions <- function(ibd_loci, phenotype, kb=500) {
  df <- ibd_loci %>%
    # Keep only relevant phenotype rows
    filter(
      if (phenotype %in% c('UC','CD')) Phenotype==phenotype | Phenotype=="IBD"
      else Phenotype %in% c("UC","CD","IBD")
    ) %>%
    # Retain necessary columns, one row per distinct signal
    dplyr::select(chr, pos, loci_ID, Signal, Phenotype) %>%
    distinct() %>%
    # Expand ± kb and sort
    mutate(start = pmax(pos - kb*1e3, 0),
           end   = pos + kb*1e3) %>%
    arrange(chr, pos) 
  
  # Convert to data.frame for an easy for-loop
  d <- as.data.frame(df)
  
  # Walk through each consecutive row (already sorted by chr, pos)
  # If expansions overlap on the same chromosome, split at midpoint
  for(i in seq_len(nrow(d) - 1)) {
    if(d$chr[i] == d$chr[i+1] && d$start[i+1] <= d$end[i]) {
      mid       <- floor((d$pos[i] + d$pos[i+1]) / 2)
      d$end[i]  <- mid
      d$start[i+1] <- mid + 1
    }
  }
  d %>% mutate(size = end-start)
}

##############
# Main
###############

if (MAKE.REGIONS){
  # Read in the GWAS regions file
  ibd.gwas.loci <- read_csv(paste0(data.dir, 'known_signals_regions_liu_delange_huang_321.csv')) %>% 
    separate(MarkerName, c("chr","pos","ref","alt"), ":",remove = FALSE) %>%
    mutate(chr = as.integer(sub("chr","",chr)), pos=as.integer(pos)) %>% 
    dplyr::select(MarkerName,chr,pos,loci_ID,Signal,Phenotype)

  # Make regions to map UC, CD and IBD GWAS variants to
  for(disease in c("UC","CD","IBD")) {
    out  <- create_non_overlapping_regions(ibd.gwas.loci, disease, kb=500)
    write.table(out, paste0(data.dir,disease,"_regions.txt"), row.names=FALSE, sep="\t", quote=FALSE)
  }
  # Make regions to make potential IBD genes to
  ibd.gene.regions <- create_non_overlapping_regions(ibd.gwas.loci, 'IBD', kb=1500)
  
  # Finding genes which fall in the IBD regions
  
  # Remove hits on X (we dont eQTL here)
  ibd.gene.regions = ibd.gene.regions[!is.na(ibd.gene.regions$chr),]
  
  # annotate this with the genes in the regions
  gene_pos = read_tsv(paste0(data.dir,"gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"))
  gene_pos$TSS = ifelse(gene_pos$strand == "+", gene_pos$start, gene_pos$end)
  ibd.gene.regions$including_genes <- sapply(1:nrow(ibd.gene.regions), function(i) {
    # Extract the current row values
    chr <- ibd.gene.regions$chr[i]
    start <- ibd.gene.regions$start[i]
    end <- ibd.gene.regions$end[i]
    # Subset genes from gene_pos based on conditions
    genes <- gene_pos[gene_pos$chromosome == chr & gene_pos$TSS > start & gene_pos$TSS < end, ]$feature_id
    # Collapse the genes into a single string
    paste(genes, collapse = ",")
  })
  ibd.gene.regions <- ibd.gene.regions %>% dplyr::select(loci_ID, Signal, including_genes)
  write_tsv(ibd.gene.regions,paste0(data.dir, 'ibd_regions.genes.txt'))

} else {
  gene.regions <- read_tsv(paste0(data.dir, 'ibd_regions.genes.txt')) %>% 
    arrange(loci_ID, Signal)
  # Initiate empty tibble
  ibd.otg.colocs <- tibble()
  for (i in 1:nrow(gene.regions)){
    row <- gene.regions[i,]
    if(is.na(row[['including_genes']])){
      next
    }
    genes <- str_split_1(row[['including_genes']],',')
    for (gene in genes){
      coloc <- colocalisationsForGene(gene) 
      if (nrow(coloc) == 0){
        next
      } else{
      coloc <- coloc %>% 
        filter(Study %in% ibd.studies,
               H4 > PP4.loose.threshold)
      if (nrow(coloc) > 0){
        print(paste(gene, 'has coloc in locus', row[['loci_ID']], row[['Signal']]))
        ibd.otg.colocs <- rbind(ibd.otg.colocs, coloc)
        }
      }
    }
  }
  ibd.otg.colocs
  
  write_tsv(ibd.otg.colocs, paste0(data.dir, 'all_otg_ibd_colocs.txt'))
  
}
