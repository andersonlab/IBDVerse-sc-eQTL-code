# Bradley September 2025
# Over-representation test of colocs with variance explained greatest in each cell-type
# module load HGI/softpack/groups/hgi/ubeR

# Packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))

# Hard code options
varex_f = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv"
gsets_gene_matrix = "response_to_reviews/gene_sets/gene_set_genes_subset.tsv.gz" # Downloaded from https://github.com/wtsi-hgi/QTLight/blob/3b0b698f5be638ae0fd7745adbed03dab18f961b/assets/data/gene_set_genes.tsv.gz
gsets_info_file = "response_to_reviews/gene_sets/gene_set_info_subset.tsv" 
sumstats.all.basedir = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'

###################
# Define fishers test function
###################
go_fish = function(pathway_genes, all_genes, sig_genes){
  pathway_genes = intersect(pathway_genes, all_genes) # Subset to the hits we find associations with
  overlap <- length(intersect(sig_genes, pathway_genes))

  contingency <- matrix(c(
    overlap,                              # in sig set & in pathway
    length(sig_genes) - overlap,          # in sig set & not in pathway
    length(pathway_genes) - overlap,      # not in sig set & in pathway
    length(all_genes) - length(sig_genes) - (length(pathway_genes) - overlap) # neither
  ), nrow = 2)

  res = fisher.test(contingency, alternative = "greater")
  resdf = data.frame(
    nsig_genes = length(sig_genes),
    nbackground_genes = length(all_genes),
    npathway_genes = length(pathway_genes),
    noverlap_genes = overlap,
    OR = res$estimate,
    pvalue = res$p.value
  )
  rownames(resdf) = NULL

  return(resdf)
}


###################
# Load
###################
varex = read.delim(varex_f) %>% 
    mutate(gene_symbol = mapIds(org.Hs.eg.db, keys = phenotype_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")) %>%
    drop_na(gene_symbol)

# Check the mapping of ense to IDs is unique:
test = varex %>% distinct(phenotype_id, gene_symbol) %>% group_by(gene_symbol) %>% summarise(count=n())
if(sum(test$count>1) == 0){
  print("All ensemble IDs are uniquely mapped")
}

# Add annotationm mapping
level.mapping <- c('0'='All Cells', 
                   '1'='Major population',
                   '2'='Cell type')

tissue.mapping <- c("ct" = "Cross-site",
                    "ti" = "Terminal ileum" ,
                    "r" = "Rectum",
                    "blood" = "Blood")

annot.mapping <- read_csv(paste0('data/all_IBDverse_annotation_mastersheet.csv'), n_max = 103) %>% 
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

varex = varex %>% left_join(annot.mapping, by=c('annotation' = 'label_machine'))

#####################
# Prepping of the gene sets
#####################
# Read in the gene sets
gene_set_genes <- read.csv(gsets_gene_matrix,
                           sep='\t',
                           header=T,
                           row.names='gene')

gene_set_info <- read.csv(gsets_info_file,
                          sep='\t',
                          header=F)
colnames(gene_set_info) = c("gene_set", "category")

rownames(gene_set_info) <- gene_set_info$gene_set ## Index by gset for later
annot_data <- gene_set_genes

# Format annot data
annot_data <- sapply(colnames(annot_data), function(x) {
  return(list(
    rownames(annot_data)[which(annot_data[[x]] == 1)]
  ))
})

#########################
# Prep the annotation with the maximum varex for each disease effector_gene, subset to want celltypes
#########################
cell_type_test = c("Blood pDCs", "cDC1", "cDC2", "Colonocyte RBFOX1+", 
                  "Colon Progenitor OLMF4+", "Colonocyte Progenitor OLMF4++Ki67+", 
                  "Colonocyte BEST4+", "Colonocyte CAECAM7++ KRT20++ IFI27++", 
                  "Epithelial Progenitor OLMF4++", "Colon Progenitor OLMF4+Ki67+", 
                  "Colonocyte CAECAM7+ KRT20+ IFI27+", "Enterocyte stem cell OLFM4+ LGR5+",
                  "Colonocyte stem cell  LGR5+", "Enterocyte stem cell MKI67+", "Enterocytes progenitor crypt OLFM4++") # Will test for significant enrichment at the level of cell-types

varex$celltype_test = ifelse(varex$label_new %in% cell_type_test, varex$label_new, NA) # Add this mapping
varex$groupedcelltype_test = ifelse(varex$label_new %in% c("cDC1", "cDC2"), "cDC_grouped", NA) # Group the DCs
varex$category_test = ifelse(varex$category_new %in% c("Colonocyte", "Stem", "Myeloid"), varex$category_new, NA) # Also test at category level


varex_max = varex %>% 
  group_by(gene_symbol) %>% 
  slice_max(rsquared)

#####################
# Manually perform over-representation test overall
#####################
# Load the maximum number of genes from the unannot all sites
expr_file = paste0(sumstats.all.basedir, "dMean__unannotated_ct_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed")
all_genes <- read_tsv(expr_file,
             col_types = cols_only(
               chr = col_character(),
               start = col_character(),
               end = col_character(),
               gene_id = col_character()
             )) %>% 
              mutate(gene_symbol = mapIds(org.Hs.eg.db, keys = gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")) %>%
              pull(gene_symbol) 
names(all_genes) = NULL
all_genes = all_genes[!is.na(all_genes)] # LOTS lost

resall = do.call(rbind, lapply(annot_data, function(pathway){
  go_fish(pathway_genes = pathway, sig_genes = unique(varex$gene_symbol), all_genes = all_genes)
})) %>% 
  mutate(fdr = p.adjust(pvalue, method="BH"))

# Have a look at nominally significant associations
resall %>% filter(pvalue < 0.05)

#####################
# Manually perform over-representation test within celltypes
#####################
col_tests = c("celltype_test", "groupedcelltype_test", "category_test")
reslist = list()
for(col in col_tests){
  level = gsub("_test", "", col)
  print(paste0("..Testing at the '", level, "' level"))

  # Subset varex for this column
  varex_col = varex[!is.na(varex[,col]),]  

  # Get the different levels
  annotations = unique(varex_col[,col])
  print(paste0("..There are ", length(annotations), " annotations in this column"))


  for(annot in annotations){
    print(paste0("..Working on ", annot))

    # Get the genes with the greatest variance explained in this cell-type
    sig_genes = varex_max %>% 
      filter(!!sym(col) == !!annot) %>% 
      pull(gene_symbol) %>% 
      unique()
    
    names(sig_genes) = NULL

    # Loop through pathways
    res_annot = do.call(rbind, lapply(annot_data, function(pathway){
      go_fish(pathway_genes = pathway, sig_genes = sig_genes, all_genes = unique(varex$gene_symbol))
    })) %>% 
      mutate(
        fdr = p.adjust(pvalue, method="BH"),
        group_type = col,
        annotation = annot
      ) %>%
      rownames_to_column("pathway") %>% 
      dplyr::select(group_type, annotation, pathway, nsig_genes, nbackground_genes, npathway_genes, noverlap_genes, OR, pvalue, fdr)
    
    if(length(reslist) == 0){
        reslist[[1]] = as.data.frame(res_annot)
      } else {
        reslist[[length(reslist) +1]] = as.data.frame(res_annot)
    }
  }
}

res = do.call(rbind, reslist)

# Save
write.csv(res, "eqtl_out/coloc_GO_FISH_notch_wnt.csv", row.names=F, quote=F)

####################
# Use cluster profiler
####################
library(clusterProfiler)
library(org.Hs.eg.db)  #
col_tests = c("celltype_test", "groupedcelltype_test", "category_test")
cplist = list()

for(col in col_tests){
  level = gsub("_test", "", col)
  print(paste0("..Testing at the '", level, "' level"))

  # Subset varex for this column
  varex_col = varex[!is.na(varex[,col]),]  

  # Get the different levels
  annotations = unique(varex_col[,col])
  print(paste0("..There are ", length(annotations), " annotations in this column"))


  for(annot in annotations){
    print(paste0("..Working on ", annot))

    # Get the genes with the greatest variance explained in this cell-type
    sig_genes = varex_max %>% 
      filter(!!sym(col) == !!annot) %>% 
      pull(gene_symbol) %>% 
      unique()
    
    names(sig_genes) = NULL
    
    tryCatch({
      ego <- enrichGO(gene          = sig_genes,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "SYMBOL",
                      ont           = "BP",   
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

      egores = ego@result
      egores$group_type = col
      egores$annotation = annot

      if(length(cplist) == 0){
          cplist[[1]] = as.data.frame(egores)
        } else {
          cplist[[length(cplist) +1]] = as.data.frame(egores)
      }
    }, error = function(e) {
      message("⚠️ enrichGO failed for ", annot, " (", col, "): ", e$message)
    })
  }
}
cpall = do.call(rbind, cplist)

# Re-adjust the pvalues for just notch and Wnt signalling
notch_wnt = cpall %>% 
  filter(str_detect(Description, regex("wnt|notch", ignore_case = TRUE))) %>% 
  group_by(annotation) %>% 
  mutate(specific_p.adjust = p.adjust(pvalue, method = "BH"))