# Bradley July 2025
# module load HGI/softpack/groups/hgi/ubeR
library(ggrepel)
library(patchwork)
library(tidyverse)
library(forcats)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(caret)
library(precrec)
library(psych)
#library(ggvenn)

server <- 'farm'

tissue <- 'multi_tissue'


if (server == 'farm'){
  repo.dir <- './'
  sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
  out.dir <- './eqtl_out/'
} else if (server == 'openstack'){
  repo.dir <- '/home/ubuntu/eqtl/code/TI-sc-eQTL-code/'
  sumstats.all.basedir <- paste0('/home/rstudio/eqtl/data/IBDverse/',tissue,'/pseudobulk-base/')
  out.dir <- '/home/rstudio/eqtl/plots/sc-eQTL-figs/eqtls/'
  
}

data.dir <- paste0(repo.dir,'/data/')
clump_file = paste0(data.dir, "clumped_all.txt.gz") # Clumping of leads
sumstats.interaction.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_12-multi_tissue_interaction_results/TensorQTL_eQTLS/'


# Plottings
SAVE.PLOTS <- TRUE
SAVE.FILES <- FALSE

##################
# Read in files
##################

source(paste0(repo.dir,'qtl_plot/helper_functions.R'))
source(paste0('qtl_plot/pleitropy_helpers.r'))

# Quick load eQTL
basedir = sumstats.all.basedir
base_dir = sumstats.all.basedir
sumstat.df = read_eqtls(sumstats.all.basedir, quick=T)
sig.sumstat.df <- sumstat.df %>% 
  filter(qval < 0.05) 

# Get nominal for specific gene
specific_gene = ""
if(specific_gene != ""){
  conditions = sig.sumstat.df %>% 
    filter(symbol == !!specific_gene) %>% 
    pull(annotation)
  
  chromosome = sig.sumstat.df %>% 
    filter(symbol == !!specific_gene) %>% 
    pull(variant_id) %>% 
    head(1) %>% 
    str_split_1(":") %>% 
    head(1)
  chromosome = gsub("chr", "", chromosome) 
  phenotype = sig.sumstat.df %>% 
    filter(symbol == !!specific_gene) %>% 
    pull(phenotype_id) %>% 
    unique()
  
  cellmeta = sig.sumstat.df %>% 
    distinct(annotation, tissue, label_new)

  # Get sumstats
  nominal = do.call(rbind, lapply(conditions, function(x){
    print(paste0("..Getting ", specific_gene, " from ", x))
    path = paste0(sumstats.all.basedir, "dMean__", x, "_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.", chromosome, ".tsv")
    temp = fread(path) %>% 
      filter(phenotype_id == !!phenotype) %>% 
      mutate(
        symbol = !!specific_gene,
        annotation = !!x
      ) %>% 
      left_join(cellmeta)
  }))

  # Save
  write.table(nominal, gzfile(paste0("temp/", specific_gene, "-nominal.txt.gz")), sep = "\t", row.names = FALSE, quote = FALSE)
}


# Get sig gene x condition pairs
sig.gene.condition = sig.sumstat.df %>% 
    mutate(pheno_annotation = paste0(phenotype_id, "-", annotation)) %>% 
    pull(pheno_annotation) %>% 
    unique()

# Load in conditional sumstats, subset for non-primary signals for genes with significant effect and merge with the sig sumstats
cond_all = read_conditional_eqtls(sumstats.all.basedir) %>% 
  mutate(pheno_annotation = paste0(phenotype_id, "-", annotation)) %>% 
  filter(
    pheno_annotation %in% sig.gene.condition # Take genes where primary effect was significant
  )

# Load the clumps 
clumps = read.delim(clump_file) %>% 
  mutate(
    phenotype_variant_id = paste0(phenotype_id, "-", variant_id),
    phenotype_clump_index = paste0(phenotype_id, "-", qtl_clump_index),
)


# Attach the clumps to both the qvalue results and conditional
sig.sumstat.df = sig.sumstat.df %>% 
  mutate(phenotype_variant_id = paste0(phenotype_id, "-", variant_id)) %>% 
  left_join(clumps %>% select(phenotype_variant_id, phenotype_clump_index))

cond_all = cond_all %>% 
  mutate(phenotype_variant_id = paste0(phenotype_id, "-", variant_id)) %>% 
  left_join(clumps %>% select(phenotype_variant_id, phenotype_clump_index))

# load the ieQTLs and do the same
ieqtl = read_ieqtls(sumstats.interaction.basedir) %>% 
    filter(
        pval_adj_bh < 0.05
    ) %>% 
    mutate(
        variant_id = factor(variant_id, levels=unique(variant_id)),
        interaction_new = interaction.mapping[interaction]
    )

# Get sig, map to clumps
ieqtl = ieqtl %>% 
  mutate(
    phenotype_variant_id = paste0(phenotype_id, "-", gsub("\\_", "\\:", variant_id))
  ) %>% 
  left_join(
    clumps %>% select(phenotype_variant_id, phenotype_clump_index),
    by="phenotype_variant_id",
    relationship="many-to-many"
  )

# Save a map of clump to min res, to be used in the coloc plotting
# If an ieQTL does NOT overlap with a base eQTL, then label it as 'ieQTL' rather than the min res
base_eqtls = unique(cond_all$phenotype_clump_index)
clumps_res = cond_all %>% 
  bind_rows(ieqtl %>% 
    rowwise() %>%
    mutate(
      annotation_type = ifelse(phenotype_clump_index %in% base_eqtls, annotation_type, "ieQTL")
    ) %>% 
    ungroup()
  ) %>% 
  group_by(phenotype_clump_index) %>% 
  slice_min(Level, with_ties=F) %>% 
  distinct(phenotype_clump_index, annotation_type)

write.table(clumps_res, gzfile(paste0(out.dir, "eQTL_minres_map.txt.gz")), sep = "\t", quote=F, row.names=F)

# Also save a version where we count the number of annotations per res
clumps_res_count = cond_all %>% 
  bind_rows(ieqtl %>% 
    rowwise() %>%
    mutate(
      annotation_type = ifelse(phenotype_clump_index %in% base_eqtls, annotation_type, "ieQTL")
    ) %>% 
    ungroup()
  ) %>% 
  group_by(phenotype_clump_index, annotation_type) %>% 
  summarise(
    count=n()
  )


write.table(clumps_res_count, gzfile(paste0(out.dir, "eQTL_count_map.txt.gz")), sep = "\t", quote=F, row.names=F)


# Load in the colocs too
coloc.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/'
known.coloc.df <- get_colocs(coloc.dir, known_ibd_only = TRUE)
known.coloc.df = known.coloc.df %>% mutate(phenotype_variant_id = paste0(phenotype_id, "-", qtl_lead)) %>% 
  left_join(clumps %>% mutate(phenotype_variant_id = gsub("\\:", "\\_", phenotype_variant_id)) %>% select(phenotype_variant_id, phenotype_clump_index)) %>% 
  filter(PP.H4.abf>0.75) # Filter 

int.coloc.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed/'
int.known.coloc.df = get_colocs(int.coloc.dir, known_ibd_only = TRUE)
int.known.coloc.df = int.known.coloc.df %>% mutate(phenotype_variant_id = paste0(phenotype_id, "-", qtl_lead)) %>% 
  left_join(clumps %>% mutate(phenotype_variant_id = gsub("\\:", "\\_", phenotype_variant_id)) %>% select(phenotype_variant_id, phenotype_clump_index)) %>% 
  filter(PP.H4.abf>0.75) # Filter 

known.coloc.df = known.coloc.df %>% 
  bind_rows(int.known.coloc.df)

#################
# Plot the distribution of effect sizes for ieQTL and eQTLs - justify the adjustment of prior for coloc
#################
# Load in nominal sumstats for chr1 from unannotated_ti base and interaction_eQTLs with ses_cd_binned
# Define condition x interactions x chromosomes with significant interactions
sigint = ieqtl %>% 
  distinct(name, interaction, chromosome)

sigint = sigint %>% 
  filter(interaction == "ses_cd_binned")

int_nom = do.call(rbind, apply(sigint, 1, FUN=function(row){
  path = paste0(sumstats.interaction.basedir, row[['name']][1], "/OPTIM_pcs/interaction_output/", row[['interaction']][1], "/cis_inter1.cis_qtl_pairs.", gsub("chr", "", row[['chromosome']][1]), ".tsv")
  print(path)

  sig_hits = ieqtl %>%
    filter(
      name == !!row[['name']][1],
      interaction == !!row[['interaction']][1],
      chromosome == !!row[['chromosome']][1]
    ) %>% 
    pull(phenotype_id) %>% 
    unique()
  
  return(fread(path) %>% 
    filter(phenotype_id %in% sig_hits))

}))

condition="dMean__unannotated_ti_all"
base_chr1_path = paste0(base_dir, condition, "/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.1.tsv")
base = fread(base_chr1_path)

base_hits = sig.sumstat.df %>% 
  filter(name == !!condition) %>% 
  pull(phenotype_id) %>% 
  unique()

base = base %>%
  filter(phenotype_id %in% base_hits)


# Plot distribution of nominal effects on chr1
nomboth = base %>% 
  filter(phenotype_id %in% base_hits) %>% 
  select(phenotype_id, variant_id, slope) %>% 
  mutate(type="eQTL") %>% 
  bind_rows(
    int_nom %>%
      rename(slope = "b_g") %>% 
      select(phenotype_id, variant_id, slope) %>%
      mutate(type = "ieQTL")
  )

# Compute sd
nomboth %>% 
  group_by(type) %>%
  summarise(
    sd = sd(slope),
    min = min(slope), 
    max = max(slope),
    range = max - min
  )

ggplot(nomboth, aes(x = slope, fill = type)) +
  geom_density(alpha = 0.5) +   # or use geom_histogram if you prefer histograms
  scale_fill_manual(values = c("ieQTL" = "#DB8CD7", "eQTL" = "#2862c0ff")) +
  labs(
    x = "beta (chr1 only)",
    y = "Density",
    fill = "eQTL type"
  ) +
  theme_classic(base_size = 14)

ggsave(paste0(out.dir, "ieqtl_all_vs_eqtl_unannot_ti_nom_beta_chr1.png"), width = 5.5, height = 3)

# have a look at significant eQTl effects only
sig_both = sig.sumstat.df %>% 
  mutate(type = "eQTL") %>% 
  bind_rows(
    ieqtl %>%
      rename(slope = "b_gi") %>% 
      mutate(type = "ieQTL")
  )


ggplot(sig_both, aes(x = slope, fill = type)) +
  geom_density(alpha = 0.5) +   # or use geom_histogram if you prefer histograms
  scale_fill_manual(values = c("ieQTL" = "#DB8CD7", "eQTL" = "#2862c0ff")) +
  labs(
    x = "beta",
    y = "Density",
    fill = "eQTL type"
  ) +
  theme_classic(base_size = 14)

ggsave(paste0(out.dir, "ieqtl_vs_eqtl_beta.png"), width = 5.5, height = 3)

# Have a look at unfiltered
unfilt_both = sumstat.df %>% 
  mutate(type = "eQTL") %>% 
  bind_rows(
    read_ieqtls(sumstats.interaction.basedir) %>% 
      mutate(type = "ieQTL") %>% 
      rename(slope = "b_gi")
  )

ggplot(unfilt_both, aes(x = slope, fill = type)) +
  geom_histogram(alpha = 0.5) +   # or use geom_histogram if you prefer histograms
  scale_fill_manual(values = c("ieQTL" = "#DB8CD7", "eQTL" = "#2862c0ff")) +
  labs(
    x = "beta",
    y = "Frequency",
    fill = "eQTL type"
  ) +
  theme_classic(base_size = 14)

ggsave(paste0(out.dir, "ieqtl_vs_eqtl_beta_unfilt.png"), width = 5.5, height = 3)

#########################
# Deep dive into FUBP1 and opposite directions of effects
#########################
gene = "FUBP1"
gene_ens = unique(cond_all %>% filter(symbol == gene) %>% pull(phenotype_id))

qtl = read_conditional_eqtls(sumstats.all.basedir) %>%
    mutate(phenotype_variant_id = paste0(phenotype_id, "-", variant_id)) %>% 
    left_join(
      read.delim("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/ld_clump_tqtl/IBDverse_plus_GTEx_plus_Yazar_results_thresh0.5/clumped_all.txt") %>% 
      mutate(
        phenotype_variant_id = paste0(phenotype_id, "-", variant_id),
        phenotype_clump_index = paste0(phenotype_id, "-", qtl_clump_index),
      ) %>%
      select(phenotype_variant_id, phenotype_clump_index)
    ) %>% # Merge with the clumps from the grouping with YAZAR and GTEx too
    filter(symbol == !!gene)

coloc = known.coloc.df %>% 
    filter(gwas_trait %in% ibd.traits, PP.H4.abf>0.75, symbol == !!gene) 

# Which clumps do we want to test
clump_leads = c("ENSG00000162613-chr1:77979080:C:CGGCCG", "ENSG00000162613-chr1:77461459:G:A")
clump_list = vector("list", length = length(clump_leads))
for(c in 1:length(clump_leads)){
  clump_list[[c]] <- unlist(strsplit(unique(qtl[qtl$phenotype_clump_index == clump_leads[c],]$phenotype_variant_id), "\\-"))[c(F,T)]
}
names(clump_list) = clump_leads

# Which conditions do we want to test?
# CD8+ TRM TGCR2+ Blood cells, rectal ADGRL3-expressing myofibroblasts and cDC2s from the rectum
conv = sig.sumstat.df %>% distinct(label_new, annotation) # Manually select annotations
condition_labels = c("T_2_blood", "Mesenchymal_2_r", "Myeloid_8_r")

# Load the genotypes
# Filter this in another window to a tempdir (much faster)
write.table(unlist(strsplit(unique(qtl$phenotype_variant_id), "\\-"))[c(F,T)], "temp/want_variants.txt", sep = "\n", quote=F, row.names=F, col.names=F)
# Add some variants on to the end of this (to ensure we keep the ones from other studies)
# echo -e "chr1:77979080:C:CGGCCG\nchr1:78157942:C:T" >> "temp/want_variants.txt"
# module load HGI/softpack/groups/macromapsqtl/macromapsqtlR4/10
# bcftools view -i ID=@temp/want_variants.txt /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/imputed.vcf.gz -o temp/want.vcf -O v

vcf.path = "temp/want.vcf"
vcf.df <- read.table(vcf.path,skip=50, header=TRUE, comment.char = "")
vcf.filt = vcf.df %>% 
  select(-c("X.CHROM", "POS", "REF", "ALT" ,"QUAL" ,"FILTER", "INFO", "FORMAT")) %>% 
  column_to_rownames("ID") %>% t()

samps = rownames(vcf.filt)

# Format
vcf.filt <- apply(vcf.filt, MARGIN=2, FUN=function(x){
  temp = sapply(strsplit(x, ":"), `[`, 1)
  temp = gsub("0\\|0", 0, temp)
  temp = gsub("0\\|1", 1, temp)
  temp = gsub("1\\|0", 1, temp)
  temp = gsub("1\\|1", 2, temp)
  return(as.numeric(temp))
}) %>% 
  as.data.frame() %>% 
  mutate(Genotyping_ID = samps)


# For each condition, load in the expression
expr = NULL
for(x in condition_labels){
  f=paste0(base_dir, "dMean__", x, "_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed")
  temp = read.delim(f) %>% 
    select(-c("X.chr", "start", "end")) %>% 
    filter(gene_id == !!gene_ens) %>% 
    column_to_rownames("gene_id") %>% 
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("Genotyping_ID")
  
  tissue=unlist(strsplit(x, "\\_")) %>% tail(1)
  format_name = paste0(
    conv[conv$annotation == x,]$label_new, " - ", tissue
  )
  format_name = make.names(gsub("\\ ", "_", format_name))

  colnames(temp)[2] = format_name
  
  if(is.null(expr)){
    expr = temp
  } else {
    expr = base::merge(expr, temp, by="Genotyping_ID", all=T)
  }
}

# Combine
both = merge(vcf.filt, expr, by="Genotyping_ID") %>% 
  column_to_rownames("Genotyping_ID")

# Manually add variants to the clump list (to include the variants from other studies)
clump_list[['ENSG00000162613-chr1:77979080:C:CGGCCG']] = c(clump_list[['ENSG00000162613-chr1:77979080:C:CGGCCG']], "chr1:77979080:C:CGGCCG", "chr1:78157942:C:T")

# Test LM
format_conds = colnames(both[,-grep("chr", colnames(both))])
res = NULL
for(c in format_conds){
  for(v in unlist(clump_list)){
      model <- lm(formula(paste0("`", c, "` ~ `", v, "`")), data = both)
      model_summary <- summary(model)
      betas <- model_summary$coefficients[, "Estimate"][-1]
      pvals <- model_summary$coefficients[, "Pr(>|t|)"][-1]
      temp = data.frame(variant = v, condition=c, gene=gene, beta=betas, pval = pvals)
      rownames(temp) = NULL
      if(is.null(res)){
        res = temp
      } else {
        res = rbind(res, temp)
      }
  }
}

# Annotate the clumps 
x <- unlist(clump_list)
names(x) <- sub("\\d+$", "", names(x))
clump_map <- data.frame(
  clump_index = names(x),
  variant = unlist(x),
  row.names = NULL
) %>% 
  rowwise() %>%
  mutate(clump_index = unlist(strsplit(clump_index, "\\-"))[c(F,T)])

colnames(clump_map)[2] = "variant"
res = merge(res, clump_map)
res = res %>% 
  mutate(
    is_sig = pval<0.05,
    is_positive = beta > 0,
    log10_pval = -log10(pval)
  )

# Re-annotate the index
res$clump_index = gsub("chr1:77979080:C:CGGCCG", "chr1:77984833:C:A", res$clump_index)

# Plot per condition
ggplot(res, aes(x = beta, y = log10_pval, fill = clump_index)) + 
  geom_point(shape = 21, color = "black", size = 3) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "lightgrey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
  facet_wrap(~ condition, nrow = 1) + 
  scale_fill_manual(values = c("chr1:77984833:C:A" = "darkorange2", "chr1:77461459:G:A" = "steelblue")) +  
  labs(
    title = paste("Nominal eQTL effects"),
    x = "Beta",
    y = expression(-log[10](p-value)),
    fill = "Clump Index"
  ) + 
  theme(legend.position="bottom") + 
  theme_bw()

ggsave(paste0(out.dir,gene, "-check_other_eQTLs.png"), width = 10, height = 4)


# Conditioning the effect of the lead eQTL on the second (chr1:77461459:G:A)
res_cond = NULL
for(c in format_conds){
  for(v in unlist(clump_list)){
    if(v != "chr1:77461459:G:A"){
      model <- lm(formula(paste0("`", c, "` ~ `", v, "` + `chr1:77461459:G:A`")), data = both)
      model_summary <- summary(model)
      betas <- model_summary$coefficients[, "Estimate"][2]
      pvals <- model_summary$coefficients[, "Pr(>|t|)"][2]
      temp = data.frame(variant = v, condition=c, gene=gene, beta=betas, pval = pvals)
      rownames(temp) = NULL
      if(is.null(res)){
        res_cond = temp
      } else {
        res_cond = rbind(res_cond, temp)
      }
    }
  }
}

res_cond = merge(res_cond, clump_map)
res_cond = res_cond %>% 
  mutate(
    is_sig = pval<0.05,
    is_positive = beta > 0,
    log10_pval = -log10(pval)
  )

# Re-annotate the index
res_cond$clump_index = gsub("chr1:77979080:C:CGGCCG", "chr1:77984833:C:A", res_cond$clump_index)

# Plot after conditioning
ggplot(res_cond, aes(x = beta, y = log10_pval, fill = clump_index)) + 
  geom_point(shape = 21, color = "black", size = 3) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "lightgrey") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
  facet_wrap(~ condition, nrow = 1) + 
  scale_fill_manual(values = c("chr1:77984833:C:A" = "darkorange2", "chr1:77461459:G:A" = "steelblue")) +  
  labs(
    title = paste("Nominal eQTL effects - condition on chr1:77461459:G:A"),
    x = "Beta",
    y = expression(-log[10](p-value)),
    fill = "Clump Index"
  ) + 
  theme(legend.position="bottom") + 
  theme_bw()

ggsave(paste0(out.dir,gene, "-check_other_eQTLs_conditional.png"), width = 10, height = 4)

#########################
# Comparing the intersection of eQTLs and ieQTLs at various LD thresholds
#########################
# Load eQTLs
ieqtl = read_ieqtls(sumstats.interaction.basedir) %>% 
    filter(
        pval_adj_bh < 0.05
    ) %>% 
    mutate(
        variant_id = factor(variant_id, levels=unique(variant_id)),
        interaction_new = interaction.mapping[interaction],
        phenotype_variant_id = paste(phenotype_id, variant_id, sep = "-")
    )

# Load clumps for each threshold
clump_dir = "../ld_clump_tqtl"
thresholds = as.character(seq(0.1, 0.5, by = 0.1))
intersection = NULL
for(x in thresholds){
  print(paste0("..Threshold: ", x))
  # Load in clumps
  temp_clumps = read.delim(paste0(clump_dir, "/results_thresh", x, "/clumped_all.txt")) %>% 
    mutate(
      phenotype_variant_id = paste(phenotype_id, variant_id, sep = "-"),
      phenotype_clump_index = paste(phenotype_id, qtl_clump_index, sep = "-")
    ) %>% 
  select(phenotype_variant_id, phenotype_clump_index)
  colnames(temp_clumps)[2] = paste0("threshold_", x)

  # Make copy of original
  temp_cond = cond_all
  temp_ieqtl = ieqtl

  # Merge onto eQTLs and ieQTLs
  temp_ieqtl = temp_ieqtl %>% left_join(
    temp_clumps, by="phenotype_variant_id", relationship="many-to-many")

  temp_cond = temp_cond %>% 
    left_join(temp_clumps, by="phenotype_variant_id", relationship="many-to-many")

  # ieQTL clumps 
  ieqtls = unique(temp_ieqtl[[paste0("threshold_", x)]])

  # eQTL clumps
  eqtls = unique(temp_cond[[paste0("threshold_", x)]])

  # ieQTL that are eQTL
  shared = intersect(ieqtls, eqtls)

  # Produce result df
  tempres = data.frame(
    threshold = x,
    number_ieQTLs = length(ieqtls),
    shared_eQTLs = length(shared)
  )

  if(is.null(intersection)){
    intersection=tempres
  } else {
    intersection = rbind(intersection, tempres)
  }

}
intersection$proportion_shared = 100*intersection$shared_eQTLs/intersection$number_ieQTLs

# Plot
ggplot(intersection, aes(x = threshold, y = proportion_shared)) + 
  geom_bar(stat = "identity", fill = "lightgrey") +  
  geom_text(aes(label = paste0(round(proportion_shared, 1), "%")), 
            vjust = -0.5, size = 3.5) +
  ylim(0, 100) + 
  labs(
    x = expression("LD clumping threshold" ~(r^2), ")"),
    y = "Intersection of ieQTLs with eQTLs (%)"
  ) + 
  theme_classic()

ggsave(paste0(out.dir,"ieqtl_vs_eqtl_intersect_threshold.png"), width = 6, height = 4)


######################
# Comparing disease effector genes from gut enriched immune cells with those of Cuomo et al., 2025
######################
# Make sure we do PP.H4> 0.8, gut and immune cell only from ours
# Cuomo et al., (2025) 
cuomo_res = read.csv("temp/cuomo_st5.csv") %>% 
  filter(
    PP.H4. > 0.8, 
    Trait.name %in% c("Inflammatory bowel disease")
  )   

cuomo_genes <- cuomo_res %>% 
  pull(Gene..Ensembl.ID.) %>% 
  unique()

# IBDverse - gut only 
ibdverse_match = known.coloc.df %>% 
  filter(
    PP.H4.abf > 0.8, 
    tissue %in% c("Rectum", "Terminal ileum"),
    category_new %in% c("B", "Plasma", "Myeloid", "T"),
    gwas_trait == "IBD"
  ) 
  
ibdverse_match_genes = ibdverse_match %>% 
  pull(phenotype_id) %>% unique()

# Count the number of genes in IBDverse
length(ibdverse_match_genes) # 78
# Count how many also nominated by Cuomo et al
length(intersect(ibdverse_match_genes, cuomo_genes)) # 38
# Proportion found previously?
length(intersect(ibdverse_match_genes, cuomo_genes))/length(ibdverse_match_genes)



# Map to loci based on genes (treat the TSS as gwas lead pos)
pheno_pos = read.delim(paste0(data.dir,  "gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt")) %>% 
  mutate(
    Gene..Ensembl.ID. = feature_id,
    TSS = ifelse(strand == "+", start, end),
    ) %>% 
  select(-c("feature_id"))

cuomo_res = cuomo_res %>% 
  left_join(pheno_pos) %>% 
  mutate(
    chr=as.numeric(chromosome),
    gwas_lead_pos=TSS,
    region_file="IBD"
  ) 

cuomo_res = map_var_to_region(cuomo_res)

# Compare the loci nominated by each
cuomo_loci = unique(cuomo_res$loci_ID)
ibdverse_loci = unique(ibdverse_match$loci_ID)
# novel
length(setdiff(ibdverse_loci, cuomo_loci))

# Plot as pie
loci.321.pi.df <- tibble(loci_ID = cuomo_loci, source='Cuomo',
                         label = 'Cuomo et al.,\n(2025)') %>% 
  bind_rows(tibble(loci_ID = ibdverse_loci, source='IBDverse',
                   label = '\n\n  GI-immune\n  IBDverse')) %>%
  bind_rows(tibble(loci_ID = ibd.regions$loci_ID, source='No colocs',
                   label = 'Remaining loci')) %>%
  distinct(loci_ID,.keep_all = TRUE) %>% 
  add_count(source) %>% 
  mutate(source = factor(source, levels = c('No colocs','IBDverse','Cuomo')))#,
         # label = paste0(label,'\nn=',n))

ggplot(loci.321.pi.df, aes(x = '', fill = source)) +
  geom_bar(colour='black') +
  geom_text(stat = 'count',aes(label = paste0(label,'\nn=',n),
                               y=after_stat(count),colour=source),
            position = position_stack(vjust = 0.5),
            size=4.3,fontface='bold') +
  coord_polar("y", start = 0) + 
  # ylim(0, 321) +
  scale_y_continuous(breaks = cumsum(unique(loci.321.pi.df$n))) +
  scale_fill_manual(values=c('IBDverse' = "#F56289",
                             'Cuomo' =  "#F099AF",
                             'No colocs' = 'white')) +
  scale_color_manual(values=c('IBDverse' = 'black',
                             'Cuomo' = 'black',
                             'No colocs' = 'black')) +
  guides(fill='none', colour='none') +
  theme_void() +
  theme(axis.text = element_text())# +
  # theme(axis.text = element_blank())

ggsave(paste0(out.dir,"IBDverse_Cuomo.png"), width = 8, height = 4)


######################
# Comparing the disease effector genes we get compared to otar. Can we find anything that drives this?
######################
## Gene per gene comparison
# Load OTAR coloc genes
otar.ibd.coloc.path <- paste0(repo.dir, 'data/all_otg_ibd_colocs.txt')
otar.ibd.colocs <- read_tsv(otar.ibd.coloc.path) %>% 
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

otar_coloc_genes <- otar.ibd.colocs %>% 
  pull(Gene_symbol) %>% 
  unique()

# How many IBDverse coloc genes are novel
length(unique(known.coloc.df$phenotype_id))
ibdverse_only = setdiff(unique(known.coloc.df$phenotype_id), unique(otar_coloc_genes))
length(ibdverse_only) # 257 

# Add Cuomo genes 
cuomo_genes = read.csv("temp/cuomo_st5.csv") %>% 
  filter(
    PP.H4. > 0.8, 
    Trait.name %in% c("Inflammatory bowel disease")
  ) %>% 
  pull(Gene..Ensembl.ID.) %>% 
  unique()

benchmark = unique(c(otar_coloc_genes, cuomo_genes))

# How many genes in IBDverse, how many novel
total_ibdverse = length(unique(known.coloc.df$phenotype_id))
novel_ibdverse = setdiff(known.coloc.df$phenotype_id, benchmark)

# How many genes are detected in blood populations of IBDverse only?
blood_genes = known.coloc.df %>% 
  filter(tissue == "Blood") %>% 
  pull(phenotype_id)

blood_only_genes = setdiff(blood_genes, known.coloc.df %>% filter(tissue != "Blood") %>% pull(phenotype_id) %>% unique())
length(blood_only_genes) # 3

blood_only_genes_novel = setdiff(blood_only_genes, benchmark)
known.coloc.df %>% distinct(phenotype_id, symbol, PP.H4.abf, gwas_trait, label_new) %>% group_by(phenotype_id, gwas_trait) %>% slice_min(PP.H4.abf) %>% filter(phenotype_id %in% blood_only_genes_novel)

# How many novel genes do we find in gut (keep cross-site if using checking B,Myeloid,T as these have lots of blood pops.)
gut_only_genes = known.coloc.df %>% 
  filter(tissue != "Blood") %>% 
  mutate(keep = case_when(
    tissue == "Cross-site" & category_new %in% c("B", "Myeloid", "T") ~ F,
    TRUE ~ TRUE)
  ) %>% 
  filter(keep) %>% 
  pull(phenotype_id) %>% 
  unique()


gut_novel_genes = setdiff(gut_only_genes, benchmark)

###############
# More formal analysis - detection vs disease sample proportion
###############
# For each annotation in the TI (including the all cells level, quantify the proportion of samples derived from CD samples)
ti_annots = known.coloc.df %>% 
  filter(tissue == "Terminal ileum") %>% 
  pull(condition_name) %>% 
  unique()

# Remove the interactions for this analysis
ti_annots = ti_annots[-grep("_cd_|_sex|_disease_status", ti_annots)]

# disease_status metadata
meta = read.delim("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/interaction_files/2024_12_27-multi_tissue/eQTL_interactions.tsv")

# Quantify the proportion of CD samples per annotation
# Quantify the number of eGenes and disease effector genes
metadir = paste0(sumstats.all.basedir, "../metadata/")
allres = do.call(rbind, lapply(ti_annots, function(x){
  print(x)
  annotmeta = read.delim(paste0(metadir, gsub("_all", "", x), "___genotype_phenotype_mapping.tsv")) %>% 
    rename(Genotyping_ID = Genotype ) %>% 
    left_join(meta, by="Genotyping_ID")

  # Get prop CD
  prop_cd = sum(annotmeta$disease_status == 1) / nrow(annotmeta)
  totaln = nrow(annotmeta) 

  # Get n eGenes
  egenes = sig.sumstat.df %>% 
    filter(name == !!x) %>% 
    pull(phenotype_id) %>% 
    unique()
  negenes = length(egenes)

  # Get colocs
  n_colocs = known.coloc.df %>% 
    filter(condition_name == !!x) %>% 
    pull(phenotype_id) %>% 
    unique() %>% 
    length()  

  # Get nice name
  label_new = known.coloc.df %>% 
    filter(condition_name == !!x) %>% 
    pull(label_new) %>% 
    unique()

  res = data.frame(condition_name = x, label_new = label_new, total = totaln, prop_cd = prop_cd, negene = negenes, nde = n_colocs, or = n_colocs / negenes)
  return(res)
}))

# Get the tissue level out to normalise
tissue_prop = allres[grep("unannotated", allres$condition_name),]$prop_cd
allres = allres %>% 
  rowwise() %>% 
  mutate(cd_enr = prop_cd / tissue_prop)

# Merge with the annotations to colour points
allres = allres %>% 
  left_join(
    known.coloc.df %>% 
      distinct(condition_name, annotation_type, category_new)
  ) %>% 
  mutate(annotation_type = factor(annotation_type, levels = c("All Cells", "Major population", "Cell type")))

# Fit linear model
lm_fit <- lm(log10(or) ~ cd_enr, data = allres)
slope <- signif(coef(lm_fit)[["cd_enr"]], 3)
pval <- signif(summary(lm_fit)$coefficients["cd_enr", "Pr(>|t|)"], 3)
# Print
print(paste0("Slope: ", slope, ". p=", pval))

# Get names of those with top5 enrichment
top5 = allres %>% 
  arrange(-or) %>% 
  pull(label_new) %>% 
  head(5)

allres = allres %>%
  rowwise() %>% 
  mutate(label = ifelse(label_new %in% top5, label_new, ""))

# Plot
ggplot(allres, aes(x = cd_enr, y = log10(or))) + 
  geom_point(aes(fill = category_new, size = annotation_type), shape = 21, stroke = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
  scale_size_manual(values = annot.class.sizes, name = 'Annotation granularity') +
  scale_fill_manual(values = umap.category.palette, name = 'Major population') +
  labs(
    title = "Disease effector gene detection vs CD sample proportion (TI)",
    x = "Enrichment of CD samples per annotation\n(normalised to 'All Cells')",
    y = expression(log[10]("Odds ratio")),
    fill = "Major population"
  ) +
  guides(size = guide_legend(order = 1), 
         fill = guide_legend(override.aes = list(size = 5), order = 0, ncol = 2)) + 
  theme_classic()

ggsave(paste0(out.dir,"CD_enrichment_vs_DE-OR.png"), width = 8, height = 4)


ggplot(allres, aes(x = prop_cd, y = log10(or))) + 
  geom_point(aes(fill = category_new, size = annotation_type), shape = 21, stroke = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
  scale_size_manual(values = annot.class.sizes, name = 'Annotation granularity') +
  scale_fill_manual(values = umap.category.palette, name = 'Major population') +
  labs(
    title = "Disease effector gene detection vs CD sample proportion (TI)",
    x = "Proportion of CD samples per annotation",
    y = expression(log[10]("Odds ratio")),
    fill = "Major population"
  ) +
  guides(size = guide_legend(order = 1), 
         fill = guide_legend(override.aes = list(size = 5), order = 0, ncol = 2)) + 
  theme_classic()

ggsave(paste0(out.dir,"CD_enrichment_vs_DE-OR_unnormalised.png"), width = 8, height = 4)


######################
# Compare the depth and nCells for cells from the TI across disease status
#####################
dat = read.csv(paste0(out.dir,"per_sample_data_ti_atlas.csv"), row.names=1)
df_long <- dat %>%
  pivot_longer(
    cols = c(Median_nCounts, ncells),  # columns to pivot
    names_to = "variable",                # name for the variable column
    values_to = "value"                   # name for the values column
  ) %>% 
  mutate(
    disease_status = factor(disease_status, levels = c("Healthy", "CD")),
    variable = factor(variable, levels = c("ncells", "Median_nCounts"))
  )


# Define custom colors
my_colors <- c("Healthy" = "navy", "CD" = "darkgoldenrod2")

# Create the plots
for(var in unique(df_long$variable)){
  ylab = ifelse(var == "ncells", "Cells per sample", "Median total counts per sample")
  ggplot(df_long %>% filter(variable == !!var), aes(x = disease_status, y = value, fill = disease_status)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
    stat_compare_means(method = "wilcox.test", 
                      comparisons = list(c("Healthy", "CD")),
                      label = "p.format") +
    scale_fill_manual(values = my_colors) +
    labs(x = "Sample type", y = ylab, fill = "") + 
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 12),
          strip.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = "none"
          )

    ggsave(paste0(out.dir,var,"_TI.png"), width = 3, height = 4)      
  
}


ggplot(df_long, aes(x = disease_status, y = value, fill = disease_status)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("Healthy", "CD")),
                     label = "p.format") +
  facet_grid(. ~ variable, 
             scales = "free_y",
             labeller = as_labeller(c("ncells" = "Cells per sample",
                                     "Median_nCounts" = "Median number of cells per sample"))) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Sample type", y = "", fill = "") + 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

ggsave(paste0(out.dir,"depth_ncells_TI.png"), width = 8, height = 4)

######################
# Quantifying the number of disease effector genes we find in immune cells from gut rather than blood
######################
immune_pops = c("B", "Plasma", "Myeloid", "T")
gut = known.coloc.df %>% 
  filter(
    tissue %in% c("Terminal ileum", "Rectum"), 
    category_new %in% immune_pops
  ) %>% 
  pull(phenotype_id) %>% 
  unique()

blood = known.coloc.df %>% 
  filter(
    tissue == "Blood", 
    category_new %in% immune_pops
  ) %>% 
  pull(phenotype_id) %>% 
  unique()

gut_immune_only = setdiff(gut, blood)
length(gut_immune_only) # 80
length(gut) # 133

##
# Of those not in gut, do we get a nominal association?
varex = read.delim("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv")

# Get list of associations to test
assocs_to_test = known.coloc.df %>% 
  filter(
    tissue %in% c("Terminal ileum", "Rectum"), 
    category_new %in% immune_pops,
    phenotype_id %in% gut_immune_only
  ) %>% 
  distinct(annotation, gwas_lead, phenotype_id) %>% 
  mutate(info = paste(gwas_lead, phenotype_id, sep = "-")) %>% 
  pull(info) %>% 
  unique()

varex_blood = varex %>%
  filter(grepl("blood", annotation)) %>% 
  mutate(
    gwas_lead = gsub("\\:", "\\_", variant_id),
    info = paste(gwas_lead, phenotype_id, sep = "-")
    ) %>% 
  filter(info %in% assocs_to_test)

# How many of these genes have a nominal significant pvalue?
varex_blood %>% 
  filter(pvalue < 0.05) %>% 
  pull(phenotype_id) %>%
  unique() %>% 
  length() # 7



#################################
# Is disease effector gene expression correlated with eQTL strength?
#################################
cor.res = NULL
for (gene in unique(varex$phenotype_id)) {
  temp = varex %>%
    filter(
      phenotype_id == !!gene,
      pvalue < 0.05
    )
  
  # try-catch to handle failures gracefully
  tryCatch({
    res = cor.test(temp$median_expression, temp$rsquared)
    df = data.frame(phenotype_id = gene,
                    cor = res$estimate,
                    p = res$p.value,
                    nobs = nrow(temp))
    rownames(df) = NULL
    
    if (is.null(cor.res)) {
      cor.res = df
    } else {
      cor.res = rbind(cor.res, df)
    }
  },
  error = function(e) {
    message("Failed for ", gene,
            " | dim(temp) = ", paste(dim(temp), collapse = "x"),
            " | Error: ", e$message)
  })
} # Fails for two genes as it has too few observations (ENSG00000163467, ENSG00000080986)

# How many are correlated? (could be in either direction as eQTL could be up or downreg.)
cor.res %>% 
  filter(
    p < 0.05
  ) %>% 
  nrow() # 103

# How many in total?
nrow(cor.res) # 378

####################
# Plot a heatmap of the coloc enrichment results
####################
res = read.csv("eqtl_out/coloc_fgsea_notch_wnt.csv")

res <- as_tibble(res)

# create a flag for significant results
res <- res %>%
  mutate(sig_label = ifelse(pvalue < 0.05, "*", ""))

anysig = res %>% filter(sig_label == "*") %>% pull(annot_id) %>% unique()
res = res %>% filter(annot_id %in% anysig)

# Subset to reactome
res = res %>% filter(grepl("REACTOME", annot_id))

# Order 
np = unique(res$annot_id)[grep("NOTCH", unique(res$annot_id))]
wnt = unique(res$annot_id)[grep("WNT", unique(res$annot_id))]

# Add major population annot
res$annotation = ifelse(res$annotation %in% c("Colonocyte", "Stem", "Myeloid"), paste0("Major population: ", res$annotation), res$annotation)
myeloid_types = c("Blood pDCs", "cDC1", "cDC2", "cDC_grouped", "Major population: Myeloid")
nonmyeloid_types = res %>% 
  filter(!(annotation %in% myeloid_types)) %>% 
  pull(annotation) %>% 
  unique()

res = res %>% 
  mutate(
    annot_id = factor(annot_id, levels = c(np, wnt)),
    annotation = factor(annotation, levels = c(myeloid_types, nonmyeloid_types))
  )


vlines <- c(length(myeloid_types) + 0.5)  # line between myeloid and non-myeloid

# figure out where to draw horizontal lines (y-axis group boundaries)
# e.g. between NOTCH and WNT groups
hlines <- c(length(np[np %in% res$annot_id]) + 0.5)

# plot heatmap
ggplot(res, aes(x = annotation, y = annot_id, fill = annot_coef)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 5) +
  scale_fill_gradient2(low = "white", high = "red") +
  # add lines for group delineation
  geom_vline(xintercept = vlines, color = "black") +
  geom_hline(yintercept = hlines, color = "black") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  ) +
  labs(
    x = "Cellular annotation",
    y = "Pathway",
    fill = "Enrichment\ncoefficient"
  )

ggsave(paste0(out.dir,"coloc_enrichment_plot.png"), width = 10, height = 4)

####################
# Plot a heatmap of the coloc GO FISHER results
####################
res = read.csv("eqtl_out/coloc_GO_FISH_notch_wnt.csv")

res <- as_tibble(res) %>%
  filter(npathway_genes > 1)

# create a flag for significant results
res <- res %>%
  mutate(sig_label = ifelse(pvalue < 0.05, "*", ""))

anysig = res %>% filter(sig_label == "*") %>% pull(pathway) %>% unique()
res = res %>% filter(pathway %in% anysig)


np = unique(res$pathway)[grep("NOTCH", unique(res$pathway))]
wnt = unique(res$pathway)[grep("WNT", unique(res$pathway))]

# Add major population annot
res$annotation = ifelse(res$annotation %in% c("Colonocyte", "Stem", "Myeloid"), paste0("Major population: ", res$annotation), res$annotation)
myeloid_types = c("Blood pDCs", "cDC1", "cDC2", "cDC_grouped", "Major population: Myeloid")
nonmyeloid_types = res %>% 
  filter(!(annotation %in% myeloid_types)) %>% 
  pull(annotation) %>% 
  unique()

res = res %>% 
  mutate(
    pathway = factor(pathway, levels = c(np, wnt)),
    annotation = factor(annotation, levels = c(myeloid_types, nonmyeloid_types))
  )


vlines <- c(length(myeloid_types) + 0.5)  # line between myeloid and non-myeloid

# figure out where to draw horizontal lines (y-axis group boundaries)
# e.g. between NOTCH and WNT groups
hlines <- c(length(np[np %in% res$pathway]) + 0.5)

# plot heatmap
ggplot(res, aes(x = annotation, y = pathway, fill = -log10(pvalue))) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 5) +
  scale_fill_gradient2(low = "white", high = "red") +
  # add lines for group delineation
  geom_vline(xintercept = vlines, color = "black") +
  geom_hline(yintercept = hlines, color = "black") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  ) +
  labs(
    x = "Cellular annotation",
    y = "Pathway",
    fill = expression(-log[10](p-value))
  )

ggsave(paste0(out.dir,"coloc_GO_plot.png"), width = 10, height = 4)


##################
# Investigating the novelty of Notch/Wnt associations
##################
# Comment asked whether previous work has indicated Notch/Wnt regulators
# Load OTAR hits
otar.ibd.coloc.path <- paste0(repo.dir, 'data/all_otg_ibd_colocs.txt')
otar.ibd.colocs <- read_tsv(otar.ibd.coloc.path) %>% 
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

# Load Cuomo et al
cuomo_res = read.csv("temp/cuomo_st5.csv") %>% 
  filter(
    PP.H4. > 0.8, 
    Trait.name %in% c("Inflammatory bowel disease")
  )   

# Found before
prev = unique(cuomo_res$Gene.name, otar.ibd.colocs$Molecular_trait)

# Load the REACTOME Wnt / Notch pathway genes
gsets_gene_matrix = "response_to_reviews/gene_sets/gene_set_genes_subset.tsv.gz"
gsets_info_file = "response_to_reviews/gene_sets/gene_set_info_subset.tsv" 
gene_set_genes <- read.csv(gsets_gene_matrix,
                           sep='\t',
                           header=T,
                           row.names='gene')

want_cols = grepl("REACTOME", colnames(gene_set_genes))
gene_set_genes = gene_set_genes[,want_cols]

# Print a list of genes in Notch and Wnt signalling pathways from previous work
pathway_colocs = NULL
for(pathway in c("NOTCH", "WNT")){
  temp = gene_set_genes[,grep(pathway, colnames(gene_set_genes))]
  temp = temp[rowSums(temp)>1,]
  pathway_df = data.frame(
    pathway = pathway,
    symbol = rownames(temp)
  )
  
  otar_hits = otar.ibd.colocs %>% 
    filter(Molecular_trait %in% rownames(temp)) %>% 
    group_by(Molecular_trait) %>% 
    mutate(
      traits_format = gsub(" [EA]", "", Trait_reported),
      traits_format = case_when(
        traits_format == "Inflammatory bowel disease" ~ 'IBD',
        traits_format == "Crohn's disease" ~ 'CD',
        TRUE ~ 'UC'
      )
    ) %>% 
    arrange(traits_format) %>% 
    mutate(
      coloc_traits_OTAR = paste(unique(traits_format), collapse=",")
    ) %>% 
    slice_max(H4, n=3) %>%
    distinct(Molecular_trait, coloc_traits_OTAR, Tissue, Source) %>%
    mutate(otar_tissue_source = paste0(Tissue, "-", Source)) %>%
    rename(symbol = "Molecular_trait") %>%
    distinct(symbol, coloc_traits_OTAR, otar_tissue_source) %>% 
    group_by(symbol) %>% 
    mutate(otar_tissue_source = paste0(otar_tissue_source, collapse=",")) %>% 
    distinct()

  cuomo_hits = cuomo_res %>% 
    filter(Gene.name %in% rownames(temp)) %>%
    rename(symbol = "Gene.name") %>% 
    group_by(symbol) %>% 
    slice_max(PP.H4., n=3) %>% 
    mutate(cuomo_cell_type = paste0("Cuomo-", Cell.type)) %>% 
    distinct(symbol, cuomo_cell_type) %>% 
    group_by(symbol) %>% 
    mutate(cuomo_cell_type = paste0(cuomo_cell_type, collapse=",")) %>% 
    distinct() %>% 
    mutate(coloc_traits_CUOMO = "IBD")

  us_hits = known.coloc.df %>%
    filter(symbol %in% rownames(temp)) %>% 
    group_by(symbol) %>% 
    mutate(
      annotation_tissue = paste("IBDverse", label_new, tissue, sep = "-")
    ) %>%
    arrange(gwas_trait) %>% 
    mutate(
      coloc_traits_IBDverse = paste(unique(gwas_trait), collapse=",")
    ) %>% 
    slice_max(PP.H4.abf, n=3) %>% 
    mutate(
      IBDverse_annotation_tissue = paste(annotation_tissue, collapse=",")
    ) %>%
    distinct(symbol, coloc_traits_IBDverse, IBDverse_annotation_tissue)


  pathway_df = pathway_df %>% 
    left_join(otar_hits) %>%
    left_join(cuomo_hits) %>% 
    mutate(
      both = ifelse(is.na(otar_tissue_source) & !is.na(cuomo_cell_type), cuomo_cell_type, NA),
      both = ifelse(!is.na(otar_tissue_source) & is.na(cuomo_cell_type), otar_tissue_source, both),
      both = ifelse(!is.na(otar_tissue_source) & !is.na(cuomo_cell_type), paste0(otar_tissue_source, ",", cuomo_cell_type), both),
      both = ifelse(is.na(otar_tissue_source) & is.na(cuomo_cell_type), NA, both)
      ) %>%
    left_join(us_hits) %>%
    mutate(all = ifelse(!is.na(both) & !is.na(IBDverse_annotation_tissue), paste(both, IBDverse_annotation_tissue, sep = ","), NA))
    

  if(is.null(pathway_colocs)){
    pathway_colocs = pathway_df
  } else {
    pathway_colocs = rbind(pathway_colocs, pathway_df)
  }
    
}

# Quantify the number of genes for which there is previous evidence from each pathway
comparison = pathway_colocs %>%
  group_by(pathway) %>%
  summarise(
    total_genes=n(),
    total_coloc_previous = sum(!is.na(both)),
    IBDverse_only = sum((!is.na(IBDverse_annotation_tissue) & is.na(both))),
    perc_add = signif(100*IBDverse_only/total_coloc_previous)
  )

# Make a bar chart for each of these pathways to illustrate the percentage change in evidence
plot_data <- comparison %>%
  select(pathway, total_coloc_previous, IBDverse_only) %>%
  pivot_longer(cols = c(total_coloc_previous, IBDverse_only),
               names_to = "type",
               values_to = "value") %>%
  mutate(
    type = recode(type,
                       total_coloc_previous = "Open Targets / Cuomo",
                       IBDverse_only = "IBDverse only")
  ) %>% 
  group_by(pathway) %>%
  mutate(
    old = value[type == "Open Targets / Cuomo"],
    new = value[type != "Open Targets / Cuomo"],
    perc_add = paste0(100*(signif(new/old,3)), "%"),
    perc_add = ifelse(type != "Open Targets / Cuomo", perc_add, "")
  ) %>% 
  select(-c(old,new)) %>% 
  group_by(pathway) %>%
  mutate(ypos = cumsum(value) - value/2)

  

ggplot(plot_data, aes(x = pathway, y = value, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = ypos, label = perc_add),
            color = "black") +
  scale_fill_manual(values = c("Open Targets / Cuomo" = "lightgrey",
                               "IBDverse only" = otar.palette[["FALSE"]])) +
  labs(x = "Pathway", y = "Number of disease effector genes for IBD/CD/UC per pathway",
       fill = "") +
  coord_flip() + 
  theme_classic(base_size = 14)

ggsave(paste0(out.dir,"NOTCH_WNT_novelty.png"), width = 8, height = 2.5)

#################
# Are any Notch/Wnt genes monogenic
#################
monogenic_list = read.csv(paste0(data.dir, "list_monogenic_ibd_list.csv")) %>% 
  rename(symbol = "Gene") %>% 
  left_join(pathway_colocs) %>% 
  filter(!is.na(pathway))



##################
# Check variance explained of drug targets
##################
varex = read.delim("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv")
hits = c('ITGA4', 'JAK2', 'IL23R')
hits_conv = known.coloc.df %>% 
  filter(symbol %in% !!hits) %>%
  distinct(phenotype_id, symbol)

# Have a look at the directions of effect in each celltype
check_doe = function(varex, gene){
  symbol = hits_conv %>% filter(phenotype_id == !!gene) %>% pull(symbol)
  temp = varex %>% 
    filter(
      pvalue < 0.05,  
      phenotype_id == !!gene
    ) %>% 
    arrange(rsquared)
  
  # Check if signs are the same
  sign1 = sign(temp$beta[1])
  allsame = all(sign(temp$beta) == sign1)
  if(allsame){
    overall = ifelse(sign1 == -1, "DOWN", "UP")
    print(paste0(symbol, " is ", overall, "regulated in all tests"))
  } else {
    print(paste0("There are opposite effects for ", symbol))
  }
}

for(g in hits_conv$phenotype_id){
  check_doe(varex, g)
}

# [1] "JAK2 is DOWNregulated in all tests"
# [1] "There are opposite effects for ITGA4"

# *** Plot var explained ***

#################
# Check overlap of our eQTLs vs GTEx (pi1)
#################

# *** Must have also loaded the eQTLs ***
cond_all = read_conditional_eqtls(sumstats.all.basedir) %>% 
  mutate(pheno_annotation = paste0(phenotype_id, "-", annotation)) %>% 
  filter(
    pheno_annotation %in% sig.gene.condition, # Take genes where primary effect was significant
  )

# Save a list of the variants that we wish to test for each annotation
annots = unique(cond_all$annotation)
inputdir = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/pi1_pairs/input/per_annot/"
if(!file.exists(inputdir)){
  dir.create(inputdir)
}

for(annot in annots){
  gene_variant_pair = cond_all %>% 
    filter(annotation == !!annot) %>% 
    distinct(phenotype_id, variant_id)
  
  write.table(gene_variant_pair, paste0(inputdir, annot, ".txt"), sep = "\t", quote=F, row.names=F)
}

# Also save version where we group the eQTLs by their minimum resolution
# Attach clumps
clumps = read.delim(clump_file) %>% 
  mutate(
    phenotype_variant_id = paste0(phenotype_id, "-", variant_id),
    phenotype_clump_index = paste0(phenotype_id, "-", qtl_clump_index),
)

cond_all = cond_all %>% 
  mutate(phenotype_variant_id = paste0(phenotype_id, "-", variant_id)) %>% 
  left_join(clumps %>% select(phenotype_variant_id, phenotype_clump_index))

level_seperation = cond_all %>% 
  group_by(phenotype_clump_index) %>% 
  slice_min(Level) %>%
  slice_min(pval_nominal, with_ties=F) 

levels = unique(level_seperation$annotation_type)
for(level in levels){
  gene_variant_pair = level_seperation %>% 
    filter(annotation_type == !!level) %>% 
    distinct(phenotype_id, variant_id)
  
  write.table(gene_variant_pair, paste0(inputdir, gsub("\\ ", "_", level), ".txt"), sep = "\t", quote=F, row.names=F)
}

# Save a list of external files
studies = c("QTS000038", "QTS000015") # 015 is GTEx, 038 is Onek1k
eqtlcat_meta_f = "/lustre/scratch125/humgen/resources_v2/eQTL_catalogue/dataset_metadata.tsv"
eqtl_cat = read.delim(eqtlcat_meta_f)
eqtl_meta_want = eqtl_cat %>% 
      filter(
        study_id %in% studies,
        quant_method == "ge" # gene expression only
      )
options = paste0(paste0(eqtl_meta_want$study_id, "-", eqtl_meta_want$dataset_id))

write.table(options, paste0(inputdir, "../replication_datasets.txt"), sep = "\t", quote=F, row.names=F, col.name=F)

# Then we can run the pi1 comparison pipeline to explore all combinations
# See https://github.com/BradleyH017/pi1_pairs

# Load results (ignore resolution for now)
pif = "../pi1_pairs/results/pi1_all.txt"
pires = read.delim(pif) %>% 
  rename(annotation="discovery") %>% 
  filter(!(annotation %in% c("All_Cells", "Major_population", "Cell_type")))
  

# Adjust the annotation name and the study-dataset name
conv = sig.sumstat.df %>% 
  distinct(annotation, label_new, category_new, tissue, annotation_type)

eqtl_meta_want = eqtl_meta_want %>%
  mutate(replication_name = paste0(study_label, "-", sample_group))

pires = pires %>% 
  left_join(conv) %>% 
  left_join(
    eqtl_meta_want %>% 
      select(study_id, dataset_id, replication_name)
  )

# 1. Make a heatmap for the comparison between GTEx and our cross-site cel-types
pi_plt = pires %>% 
  filter(
    tissue == "Cross-site",
    annotation_type == "Cell type",
    grepl("GTEx", replication_name)
  ) %>% 
  select(pi1, replication_name, label_new) %>%
  pivot_wider(names_from = replication_name, values_from = pi1, values_fill = 0) %>%
  column_to_rownames("label_new") %>%
  as.matrix() 


# Scale and filter for complete rows
pi_plt_scaled = t(scale(t(pi_plt)))
rownames(pi_plt) <- rownames(pi_plt)
pi_plt_scaled = pi_plt_scaled[complete.cases(pi_plt_scaled),]

# Make row annotation
plot_cats = conv %>% 
    filter(
      label_new %in% rownames(pi_plt_scaled),
      tissue == "Cross-site") %>% 
    pull(category_new)

row_ha <- rowAnnotation(
  Category = plot_cats,
  col = list(Category = umap.category.palette),
  show_annotation_name = FALSE, 
  width = unit(4, "mm")         
)

ppi=300
png(paste0(out.dir, "eqtl_pi1_heatmap_gtex.png"), res=ppi, width=17*ppi, height=17*ppi)
Heatmap(pi_plt_scaled,
        name = "π1 (scaled)",
        cluster_rows = TRUE,   
        cluster_columns = TRUE, 
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_km = 2,
        column_km=2,
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"), 
        row_title = "IBDverse cell-type (cross-site)",
        column_title = "GTEx Tissue",
        left_annotation = row_ha,
        column_names_rot = 45)
dev.off()

# 2. Compare replication rates in the 'All Cells' tissue with each GTEx dataset
top_per_tissue = pires %>% 
  filter(
    tissue != "Cross-site",
    label_new == "All Cells",
    grepl("GTEx", replication_name)
  ) %>% 
  group_by(tissue) %>% 
  mutate(rank = rank(-pi1, ties.method = "first")) %>% 
  rowwise() %>%
  mutate(
    variable = as.character(replication_name),
    label = ifelse(rank <=3, variable, ""),
    label = ifelse(grepl("colon", variable), variable, label),
    label = ifelse(grepl("small_intestine", variable), variable, label),
    label = gsub("GTEx-", "", label),
    label = ifelse(label != "", paste0(label, " (", rank, ")"), label)
  ) %>%
  ungroup() %>%
  mutate(tissue = factor(tissue, levels = c("Blood", "Terminal ileum", "Rectum")))

# Plot
ggplot(top_per_tissue, aes(x=tissue, y=pi1, fill=tissue)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + 
    geom_point() + 
    scale_fill_manual(values = tissue.palette) + 
    theme_classic() + 
    geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) + 
    theme(
      legend.position="none",
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
    ) + 
    labs(
      x="'All Cells' tissue (IBDverse)",
      y="Replication rate with each GTEx tissue (π1)"
    )

ggsave(paste0(out.dir,"top_gtex_tissue_per_all_cells_tissue_pi1.png"), width = 8, height = 6)


# 3. Specific comparison of replication rates with blood and colon-transvers
resmat_spec <- pires %>%
  filter(
    replication_name %in% c('GTEx-blood', 'GTEx-colon_transverse'),
    tissue == "Cross-site",
    label_new != "All Cells",
    pi1 > 0
  ) %>%
  select(label_new, category_new, replication_name, pi1) %>% 
  mutate(
    category_new = factor(category_new, levels = c("Myeloid", "T", "B", "Plasma", "Enterocyte", "Colonocyte", "Secretory", "Stem", "Mesenchymal")),
    replication_name = gsub("GTEx-", "", replication_name)
  )
  
vals = unique(resmat_spec$replication_name)
comparisons = combn(as.character(vals), 2, simplify = FALSE)


ggplot(resmat_spec, aes(x=replication_name, y=pi1, fill=replication_name)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_classic() + 
    facet_grid(.~category_new) + 
    theme(legend.position="bottom", axis.text.x = element_text(size = 0)) + 
    stat_compare_means(comparisons = comparisons, method = "wilcox.test") + 
    labs(
      x="Major population",
      y="Replication rates (π1)",
      fill="GTEx tissue"
    )

ggsave(paste0(out.dir,"blood_vs_colon_sharing_gtex_pi1.png"), width = 10, height = 5)

# 4. Compare replication rates across resolution with GTEx
piresres = read.delim(pif) %>% 
  rename(annotation="discovery") %>% 
  filter(annotation %in% c("All_Cells", "Major_population", "Cell_type"))

# Plot these
piresresplt = piresres %>% 
  left_join(
    eqtl_meta_want %>% 
      select(study_id, dataset_id, replication_name)
  ) %>% 
  filter(grepl("GTEx", replication_name)) %>% 
  mutate(
    annotation_type = factor(gsub("\\_", " ", annotation), levels = c("All Cells", "Major population", "Cell type"))
  )

vals = unique(piresresplt$annotation_type)
comparisons = combn(as.character(vals), 2, simplify = FALSE)

ggplot(piresresplt, aes(x=annotation_type, y=pi1, fill=annotation_type)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + 
    geom_point() + 
    scale_fill_manual(values = annot.class.palette) + 
    theme_classic() + 
    stat_compare_means(comparisons = comparisons, method = "wilcox.test") + 
    theme(
      legend.position="none",
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
    ) + 
    labs(
      x="Resolution of eQTL detection (IBDverse)",
      y="Replication rate with each GTEx tissue (π1)"
    )

ggsave(paste0(out.dir,"gtex_tissue_p1_per_res.png"), width = 5, height = 6)

# Have a look
piresresplt %>% 
  group_by(annotation_type) %>% 
  slice_max(pi1) %>% 
  as.data.frame()

# Get medians
piresresplt %>% 
  group_by(annotation_type) %>% 
  mutate(median = median(pi1)) %>%
  distinct(annotation_type, median)
  

# 5. Look at general replication rates
# If we pick the maximum pi1 per annotation, what does the distirbution look like?
max_pi_plt = pires %>% 
  filter(
    grepl("GTEx", replication_name)
  ) %>% 
  group_by(annotation) %>% 
  slice_max(pi1, with_ties=F) %>% 
  filter(pi1 > 0)

ggplot(max_pi_plt, aes(x=pi1)) + 
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 50, fill = "lightgrey", color = "black") +
  geom_vline(xintercept = median(max_pi_plt$pi1), lty="dashed", color="darkorange") + 
  theme_classic() + 
  labs(x="Maximum observed π1 per annotation", y="Frequency") + 
  xlim(c(0,1))

ggsave(paste0(out.dir,"max_pi1_gtex_per_annot.png"), width = 5, height = 3)

# Look at bad replicators
max_pi_plt %>% 
  filter(pi1 < 0.25) %>% 
  head(5) %>% 
  as.data.frame()

##############
# Comparison with Yazar et al
##############

# Load metadata
studies = c("QTS000038", "QTS000015") # 015 is GTEx, 038 is Onek1k
eqtlcat_meta_f = "/lustre/scratch125/humgen/resources_v2/eQTL_catalogue/dataset_metadata.tsv"
eqtl_cat = read.delim(eqtlcat_meta_f)
eqtl_meta_want = eqtl_cat %>% 
      filter(
        study_id %in% studies,
        quant_method == "ge" # gene expression only
      )
    
conv = sig.sumstat.df %>% 
  distinct(annotation, label_new, category_new, tissue, annotation_type)

eqtl_meta_want = eqtl_meta_want %>%
  mutate(replication_name = paste0(study_label, "-", sample_group))



# Load if not already, filter for Onek1k
pif = "../pi1_pairs/results/pi1_all.txt"
pi1k = read.delim(pif) %>% 
  rename(annotation="discovery") %>% 
  filter(!(annotation %in% c("All_Cells", "Major_population", "Cell_type"))) %>% 
  left_join(conv) %>% 
  left_join(
    eqtl_meta_want %>% 
      select(study_id, dataset_id, replication_name)
  ) %>% 
  filter(
    grepl("OneK1K", replication_name)
  )

# Look at the overall replication rates, but divide by major population
max_pi1k = pi1k %>% 
  filter(
    annotation_type == "Cell type",
    tissue == "Cross-site"
  ) %>% 
  group_by(label_new) %>% 
  slice_max(pi1) %>%
  filter(pi1 > 0) %>% 
  mutate(category_new = factor(category_new, levels=c("B", "Myeloid", "T", "Plasma", "Colonocyte", "Enterocyte", "Secretory", "Stem", "Mesenchymal")))

median_by_category <- max_pi1k %>% 
  group_by(category_new) %>% 
  summarise(median_pi1 = median(pi1))

nonimmune_cats = c("Colonocyte", "Enterocyte", "Secretory", "Stem", "Mesenchymal")

ggplot(max_pi1k, aes(x=pi1, fill=category_new)) + 
  geom_rect(data = subset(max_pi1k, category_new %in% nonimmune_cats),
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.3, inherit.aes = FALSE) + 
  geom_histogram(bins = 50, color = "black") +
  facet_grid(category_new~.) + 
  scale_fill_manual(values = umap.category.palette) + 
  geom_vline(data = median_by_category, 
             aes(xintercept = median_pi1), 
             lty="dashed", color="black", size = 1) + 
  geom_text(
    data = median_by_category,
    aes(x = 0, y = 3, label = paste0(category_new, "\nmedian: ", signif(median_pi1, 2))),
    inherit.aes = FALSE, hjust = 0, vjust = 0, size=4
  ) + 
  theme_classic() + 
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
  ) + 
  labs(x="Maximum observed π1 per annotation", y="Frequency") + 
  xlim(c(0,1))

ggsave(paste0(out.dir,"max_pi1_onek1k_per_annot.png"), width = 5, height = 8)


# Have a look at best and worst
max_pi1k %>% arrange(-pi1) %>% head(5) %>% as.data.frame() # Best
max_pi1k %>% arrange(pi1) %>% head(5) %>% as.data.frame() # Worst


# 2. Compare replication rates of our cell-types with Yazar (pairwise)
# Annotate categories of OneK1K
k1_map = data.frame(
  replication_name = unique(pi1k$replication_name),
  OneK1K_category = c(
    "B", # "OneK1K-B_intermediate"   
    "B", # "OneK1K-B_memory" 
    "B", # "OneK1K-B_naive"
    "Myeloid", # "OneK1K-CD14_Mono"
    "Myeloid", # "OneK1K-CD16_Mono" 
    "T", # "OneK1K-CD4_CTL" 
    "T", # "OneK1K-CD4_Naive"
    "T", # "OneK1K-CD4_TCM"
    "T", # "OneK1K-CD4_TEM" 
    "T", # "OneK1K-CD8_Naive"
    "T", # "OneK1K-CD8_TCM"
    "T", # "OneK1K-CD8_TEM"     
    "unknown", # "OneK1K-HSPC" (Haematopoetic stem)
    "T", # "OneK1K-MAIT"
    "T", # "OneK1K-NK"
    "T", # "OneK1K-NK_CD56bright"
    "T", # OneK1K-NK_Proliferating
    "Plasma", # OneK1K-Plasmablast
    "Myeloid", # OneK1K-Platelet
    "T", # OneK1K-Treg
    "Myeloid", # OneK1K-cDC2
    "T", # "OneK1K-dnT"
    "T", # OneK1K-gdT
    "T" # "OneK1K-pDC"
  )
)

pi1k = pi1k %>% 
  left_join(k1_map)


# Get square version
pi_plt = pi1k %>% 
  filter(
    tissue == "Cross-site", 
    annotation_type == "Cell type",
    category_new %in% c("B", "Plasma", "Myeloid", "T")
  ) %>% 
  select(pi1, replication_name, label_new) %>%
  pivot_wider(names_from = replication_name, values_from = pi1, values_fill = 0) %>%
  column_to_rownames("label_new") %>%
  as.matrix() 

# Scale and filter for complete rows
pi_plt_scaled = t(scale(t(pi_plt)))
rownames(pi_plt_scaled) <- rownames(pi_plt)
pi_plt_scaled = pi_plt_scaled[complete.cases(pi_plt_scaled),]

# This time, we won't cluster, so order the cell-types based on category
categories_ordered = pi1k %>% 
  filter(
    tissue == "Cross-site", 
    annotation_type == "Cell type",
    category_new %in% c("B", "Plasma", "Myeloid", "T")
  ) %>%
  distinct(label_new, category_new) %>% 
  arrange(category_new)

onek_categories_ordered = pi1k %>% 
  distinct(replication_name, OneK1K_category) %>% 
  arrange(OneK1K_category)

pi_plt_scaled = pi_plt_scaled[match(categories_ordered$label_new, rownames(pi_plt_scaled)),]
pi_plt_scaled = pi_plt_scaled[,match(onek_categories_ordered$replication_name, colnames(pi_plt_scaled))]
pi_plt_scaled = pi_plt_scaled[complete.cases(pi_plt_scaled),]


# Define colour for 'unkown'
unkown_col = c("grey")
names(unkown_col) = "unknown"

# Make row annotation
plot_cats = categories_ordered %>% 
    filter(label_new %in% rownames(pi_plt_scaled)) %>%
    pull(category_new)

row_ha <- rowAnnotation(
  Category = plot_cats,
  col = list(Category = c(umap.category.palette, unkown_col)),
  show_annotation_name = FALSE, 
  width = unit(4, "mm")         
)

# Make column annotation
onek_cats = onek_categories_ordered %>% 
  pull(OneK1K_category)

col_ha <- columnAnnotation(
  Category = onek_cats,
  col = list(Category = c(umap.category.palette, unkown_col)),
  show_annotation_name = FALSE, 
  height = unit(4, "mm")         
)


ppi=300
png(paste0(out.dir, "eqtl_pi1_heatmap_OneK1K.png"), res=ppi, width=17*ppi, height=13*ppi)
Heatmap(pi_plt_scaled,
        name = "π1 (scaled)",
        cluster_rows = T,   
        cluster_columns = T, 
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_km = 4,
        column_km=5,
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"), 
        row_title = "IBDverse cell-type (cross-site)",
        column_title = "OneK1K cell-type",
        left_annotation = row_ha,
        top_annotation = col_ha,
        column_names_rot = 45)
dev.off()

# 3. Generate confusion matrix based on the top matching cell-type. Is it in the same category?
max_pi1k_immune = max_pi1k %>% 
  filter(
    category_new %in% c("B", "Plasma", "Myeloid", "T"),
    !(annotation %in% c("Myeloid_14_ct", "Myeloid_10_ct")) # Exclude these with very low eQTLs
  ) %>% 
  mutate(category_new = as.character(category_new)) %>% 
  left_join(k1_map)

# Create confusion matrix
confusion_matrix <- table(Predicted=max_pi1k_immune$OneK1K_category, Actual=max_pi1k_immune$category_new)
print("Confusion Matrix:")
print(confusion_matrix)

# Calculate accuracy metrics
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Overall Accuracy:", round(accuracy, 3))) 
# "Overall Accuracy: 0.763"

# 1. Chi-square test for association
chi_square_test <- chisq.test(confusion_matrix)
print("Chi-square test:")
print(chi_square_test)

#1] "Chi-square test:"
#
#        Pearson's Chi-squared test
#
#data:  confusion_matrix
#X-squared = 44.666, df = 9, p-value = 1.064e-06

# Create a visualization of the confusion matrix
confusion_df <- as.data.frame(confusion_matrix)

# Heatmap of confusion matrix
ggplot(confusion_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_classic() +
  labs(title = "IBDverse vs OneK1K major populations",
       x = "OneK1K",
       y = "IBDverse",
       fill = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(out.dir,"confusion_matrix_max_pi1_onek1k_categories.png"), width = 5, height = 5)


# Calculate per-class metrics
classes <- union(levels(factor(max_pi1k_immune$category_new)),
                 levels(factor(max_pi1k_immune$OneK1K_category)))

metrics <- data.frame(
  Class = classes,
  Precision = NA,
  Recall = NA,
  F1 = NA
)

for (i in seq_along(classes)) {
  cat <- classes[i]
  TP <- confusion_matrix[cat, cat]                      # true positives
  FP <- sum(confusion_matrix[cat, ]) - TP               # false positives
  FN <- sum(confusion_matrix[, cat]) - TP               # false negatives
  TN <- sum(confusion_matrix) - TP - FP - FN            # true negatives
  
  precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
  f1        <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                      NA, 2 * precision * recall / (precision + recall))
  
  metrics[i, c("Precision", "Recall", "F1")] <- c(precision, recall, f1)
}
metrics

#
#    Class Precision    Recall        F1
#1       B 0.6666667 0.5714286 0.6153846
#2 Myeloid 1.0000000 0.5555556 0.7142857
#3  Plasma 1.0000000 0.3333333 0.5000000
#4       T 0.7307692 1.0000000 0.8444444

# Plot F1
order = metrics %>% 
  arrange(F1) %>%
  pull(Class)

metrics = metrics %>% 
  mutate(Class = factor(Class, levels = order))

ggplot(metrics, aes(x=Class, y=F1, fill=Class)) + 
  geom_bar(stat="identity") + 
  theme_classic() + 
  scale_fill_manual(values=umap.category.palette) + 
  labs(
    x="Major population of IBDverse cell-type",
    y="Concordance of major population with top-replicating cell-type\nin Yazar et al (F1-score)"
  ) + 
  ylim(c(0,1)) + 
  theme(legend.position="none")
  
ggsave(paste0(out.dir,"F1_major_pop_us_vs_yazar.png"), width = 4, height = 5.5)


balanced_acc <- mean(metrics$Recall, na.rm = TRUE)
balanced_acc
# [1] 0.6150794

# Cohens kappa test (is agreement beyond chance)
kappa_test <- cohen.kappa(confusion_matrix)
kappa_val <- kappa_test$kappa
se_val <- sqrt(kappa_test$var.kappa)
z_val <- kappa_val / se_val
p_val <- 2 * (1 - pnorm(abs(z_val)))
print(paste0("Kappa: ", signif(kappa_val,2), ". p=", signif(p_val,2)))
# "Kappa: 0.602. p=2.1e-08"


