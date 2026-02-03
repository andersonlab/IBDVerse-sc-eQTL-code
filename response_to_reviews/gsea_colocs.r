# Bradley August 2025
# Targetted GSEA of Wnt and Notch pathways in colonocytes and cDCs respectively.
# Following https://github.com/wtsi-hgi/QTLight/blob/3b0b698f5be638ae0fd7745adbed03dab18f961b/bin/fgsea_ieQTLs.R#L171
# module load HGI/softpack/users/oa3/ieqtl_gsea/2 

# Packages
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(Matrix))

# Next two files were downloaded from  https://github.com/wtsi-hgi/QTLight/blob/3b0b698f5be638ae0fd7745adbed03dab18f961b/assets/data/
# Subset to include only Wnt and Notch pathway genes
#zgrep -i -E 'NOTCH|WNT' gene_set_info.tsv.gz > gene_set_info_subset.tsv # Get pathway names
#cut -f1 gene_set_info_subset.tsv \
#  | sed 's/\r$//' | sed 's/^[ \t]*//;s/[ \t]*$//' \
#  | awk 'length>0' \
#  | sort -u > ids.txt
#
#idxs=$(
#  zcat gene_set_genes.tsv.gz | head -1 | awk -F'\t' '
#    BEGIN {
#      # load desired names
#      while ((getline line < "ids.txt") > 0) {
#        gsub(/\r/, "", line); gsub(/^[ \t]+|[ \t]+$/, "", line);
#        if (line != "") want[line]=1
#      }
#      close("ids.txt")
#    }
#    {
#      pos[1]=1; n=1   # always keep first column
#      for (i=2; i<=NF; i++) {
#        col=$i; gsub(/\r/, "", col); gsub(/^[ \t]+|[ \t]+$/, "", col);
#        if (col in want) pos[++n]=i
#      }
#    }
#    END {
#      for (i=1; i<=n; i++) printf i==1 ? "%d" : ",%d", pos[i]
#    }'
#)
#zcat gene_set_genes.tsv.gz | cut -f"$idxs" | gzip > gene_set_genes_subset.tsv.gz

# Hard code options
varex_f = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv"
gsets_gene_matrix = "response_to_reviews/gene_sets/gene_set_genes_subset.tsv.gz" # Downloaded from https://github.com/wtsi-hgi/QTLight/blob/3b0b698f5be638ae0fd7745adbed03dab18f961b/assets/data/gene_set_genes.tsv.gz
gsets_info_file = "response_to_reviews/gene_sets/gene_set_info_subset.tsv" 

##############
# Read in the variance explained results. Add the annotation mappings, and specify the custom groupings of cells
###############
# Load
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

# Define custom cell mappings
cell_type_test = c("Blood pDCs", "cDC1", "cDC2", "Colonocyte RBFOX1+", 
                  "Colon Progenitor OLMF4+", "Colonocyte Progenitor OLMF4++Ki67+", 
                  "Colonocyte BEST4+", "Colonocyte CAECAM7++ KRT20++ IFI27++", 
                  "Epithelial Progenitor OLMF4++", "Colon Progenitor OLMF4+Ki67+", 
                  "Colonocyte CAECAM7+ KRT20+ IFI27+", "Enterocyte stem cell OLFM4+ LGR5+",
                  "Colonocyte stem cell  LGR5+", "Enterocyte stem cell MKI67+", "Enterocytes progenitor crypt OLFM4++") # Will test for significant enrichment at the level of cell-types
varex$celltype_test = ifelse(varex$label_new %in% cell_type_test, varex$label_new, NA) # Add this mapping
varex$groupedcelltype_test = ifelse(varex$label_new %in% c("cDC1", "cDC2"), "cDC_grouped", NA) # Group the DCs
varex$category_test = ifelse(varex$category_new %in% c("Colonocyte", "Stem", "Myeloid"), varex$category_new, NA) # Also test at category level

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

# Some genes may be weighted greater because the strength of the eQTL effect will be greater.
# To avoid this, scale the rsquared within gene also
varex = varex %>% 
  group_by(gene_symbol) %>% 
  mutate(rsq_scaled = scale(rsquared)) %>% 
  ungroup() %>% 
  as.data.frame()

####################
# Perform testing within the various groupings
####################
col_tests = c("celltype_test", "groupedcelltype_test", "category_test")
rank_types = c("rsquared", "rsq_scaled")
reslist = list()
for(col in col_tests){
  level = gsub("_test", "", col)
  print(paste0("..Testing at the '", level, "' level"))

  # Subset varex
  varex_col = varex[!is.na(varex[,col]),]

  # Get the different levels
  annotations = unique(varex_col[,col])
  print(paste0("..There are ", length(annotations), " annotations in this column"))

  for(annot in annotations){
    print(paste0("....Working on ", annot))

    # Subset
    varex_col_annot = varex_col[varex_col[,col] == annot,]

    for(rank_type in rank_types){

      print(paste0("......Ranking by ", rank_type))
      
      # get the maximum variance explained for each gene within the given annotation
      varex_col_annot_max = varex_col_annot %>% 
        group_by(gene_symbol) %>% 
        slice_max(!!sym(rank_type), with_ties=F)

      # Set up list
      summary_data = varex_col_annot_max[[rank_type]]
      names(summary_data) = varex_col_annot_max$gene_symbol

      pos_std = ifelse(rank_type == "rsquared", 'pos', 'std')

      # Run GSEA
      ## This will run fgsea.Multilevel, which gives more precise p-vals at a
      ## consequence for longer compute time. Could adjust to using fgsea.Simple
      ## See: https://github.com/ctlab/fgsea/blob/master/R/fgsea.R
      fgseaRes <- fgsea::fgseaMultilevel(pathways = annot_data,
                                      stats = summary_data,
                                      sampleSize = 101,
                                      eps = 0,
                                      minSize = 1,
                                      maxSize = Inf,
                                      scoreType = pos_std,
                                      nproc = 1,
                                      gseaParam = 1,
                                      BPPARAM = NULL)

      # Reformat
      fgseaRes <- fgseaRes %>%
        dplyr::rename(
          annot_id = pathway,
          annot_coef = NES,
          pvalue = pval,
          qvalue_bh = padj, ## fgsea automatically does BH correction
          n_gene_contained = size
        ) %>% 
        arrange(pvalue)

      # Add info about the test type and group
      fgseaRes$coltype = col
      fgseaRes$annotation = annot
      fgseaRes$rank_type = rank_type

      # Add to results
      if(length(reslist) == 0){
        reslist[[1]] = as.data.frame(fgseaRes)
      } else {
        reslist[[length(reslist) +1]] = as.data.frame(fgseaRes)
      }
    }
  }
}
res = do.call(rbind, reslist)

# Are there any nominally significant results? 
res %>% 
  filter(pvalue < 0.05) %>% 
  dplyr::select(annot_id, pvalue, annot_coef, annotation)

# If we just limit to the testing of REACTOME pathways and correct, what do we get? 
res %>%  
  filter(grepl("REACTOME", annot_id)) %>% 
  group_by(coltype, annotation) %>% 
  mutate(readjust = p.adjust(pvalue, method = "BH")) %>% 
  arrange(pvalue) %>%
  dplyr::select(annot_id, pvalue, readjust, annot_coef, annotation, rank_type) %>% 
  head(10)

# Save
write.csv(res %>% dplyr::select(-leadingEdge), "eqtl_out/coloc_fgsea_notch_wnt.csv", row.names=F, quote=F)