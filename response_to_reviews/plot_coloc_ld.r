# Bradley August 2025
# Plot coloc plots, colour points by LD
# module load HGI/softpack/users/bh18/Colocalisation_analysis/2
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000184384" "MAML2" "Myeloid_4_ti" "UC" "F"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000172575" "RASGRP1" "Colonocyte_ct" "CD" "F"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000145012" "LPP" "Secretory_r" "IBD" "F" # For Jingling
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000172575" "RASGRP1" "Colonocyte_ct" "CD" "T" "F" "F" # Plot the minres
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000136997" "MYC" "Epithelial_12_ct" "CD" "T" # Plot the minres
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000136997" "MYC" "Epithelial_12_ct" "CD" "F" # Don't plot the minres



################
# Preamble
################
# Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))


# Options
args = commandArgs(trailingOnly=TRUE)
gene = args[1] # gene="ENSG00000172575"
symbol = args[2] # symbol = "RASGRP1"
condition = args[3] # condition="Colonocyte_ct"
gwas_trait = args[4] # gwas_trait="CD"
plot_minres = as.logical(args[5]) # plot_minres = T # Do you want some kind of annotation for the min res of variants?
plot_minres_lines = as.logical(args[6]) # plot_minres_lines = F # Do you want the min res plotted as lines (TRUE)? Or as dots (FALSE)?
plot_eqtl_count = as.logical(args[7]) # plot_eqtl_count = F # This will plot the position and number of eQTLs per resolution on top?
print(paste0("~~~~~~ Plotting coloc between ", symbol, " and ", gwas_trait, " in ", condition))

# Hard code paths and options
sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
gwas_dir <- "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/gwas_sumstats_final/"
gwas_master_list = paste0(gwas_dir, "gwas_files.txt")
clump_file = "data/clumped_all.txt.gz"
window = 5e5 # Centred on eQTL lead, in either direction
plink_prefix = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/plink/imputed_chr" # Location of plink genotyping files
annotation_file = "data/all_IBDverse_annotation_mastersheet.csv"
out = "response_to_reviews/coloc_plots"
gene_pos_f = "data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
coloc_dir = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/'
tempout = "temp/coloc"
eqtl_minres_map = "eqtl_out/eQTL_minres_map.txt.gz"
eQTL_count_map = "eqtl_out/eQTL_count_map.txt.gz"
for(dir in c(out, tempout)){
    if(!file.exists(dir)){
        dir.create(dir)
    }
}
annot.class.palette <- c('Tissue'='#edf8b1',
                         'All Cells'='#edf8b1',
                         'Major category'='#7fcdbb',
                         'Major population'='#7fcdbb',
                         'Cell type'='#2c7fb8',
                         'ieQTL' = '#DB8CD7')

############
# Get gwas and QTL sumstats
############
## eQTL
# Lead
eqtlf = paste0(sumstats.all.basedir, "dMean__", condition, "_all/OPTIM_pcs/base_output/base/Cis_eqtls_qval.tsv")
eqtl_q = fread(eqtlf) 
eqtl_q_want = eqtl_q %>% 
    filter(phenotype_id == !!gene)
eqtl_lead = eqtl_q_want %>% 
    pull(variant_id)

split_lead = unlist(strsplit(eqtl_lead, "\\:"))
chrnum = gsub("chr", "", split_lead[1])
leadpos = as.numeric(split_lead[2])
min=max(0,leadpos-window)
max=leadpos+window

# Nominal
print("..Loading the eQTL sumstats")
eqtlnom = paste0(sumstats.all.basedir, "dMean__", condition, "_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.", chrnum, ".tsv")
eqtl = fread(eqtlnom) %>% 
    rowwise() %>% 
    mutate(pos = as.numeric(unlist(strsplit(variant_id, "\\:"))[c(F,T,F,F)])) %>% 
    filter(
        pos > min,
        pos < max,
        phenotype_id == !!gene
    )

## GWAS
print("..Loading the GWAS sumstats")
# Get long name
long_name = read.delim(gwas_master_list, header=F) %>% 
    filter(V1 == !!gwas_trait) %>%
    pull(V2)

# Load 
gwasf = paste0(gwas_dir, long_name, ".txt.gz")
gwas = fread(gwasf) %>% 
    filter(
        Chr == !!chrnum,
        Pos > min, 
        Pos < max
    )

gwas_lead = gwas %>% 
    slice_min(pval) %>% 
    pull(RSid)

############
# Compute LD for each variant with the index per study
############
# Specify plink file
plinkfile = paste0(plink_prefix, chrnum)

# Save a list of variants for each
qtlvartemp = paste0(tempout, "/", gene, "-", condition, "-", "qtlvariants.txt")
gwasvartemp = paste0(tempout, "/", gene, "-", gwas_trait, "-", "gwasvariants.txt")
write.table(eqtl$variant_id, qtlvartemp, quote=F, sep = "\t", row.names=F, col.names=F)
write.table(gwas$RSid, gwasvartemp, quote=F, sep = "\t", row.names=F, col.names=F)

# For each set, process variants
for(variantlist in c(qtlvartemp, gwasvartemp)){
    if(variantlist == qtlvartemp){
        type = "QTL"
    } else {
        type = "GWAS"
    }
    print(paste0("..Computing LD for ", type, " variants"))
    # Filter the plink file
    subset_out = paste0(tempout, "/", gene, "-", condition, "-", type)
    system(sprintf('plink --bfile %s --extract %s --make-bed --out %s --silent',  plinkfile, variantlist, subset_out))

    # Now compute LD using this
    ld_out = paste0(tempout, "/", gene, "-", condition, "-", type, "-LD")
    target = ifelse(type == "QTL", eqtl_lead, gwas_lead)
    system(sprintf('plink --bfile %s --ld-window-kb 50000000 --ld-snp %s --ld-window 20000 --ld-window-r2 0 --r2 --out %s --silent', subset_out, target, ld_out)) 
}

# Load in each 
eqtl_ld = read.delim(paste0(tempout, "/", gene, "-", condition, "-QTL-LD.ld"), sep = "") %>% 
    dplyr::rename(
        variant_id = SNP_B,
        r2_lead = R2
    )%>% 
    select(variant_id, r2_lead)

gwas_ld = read.delim(paste0(tempout, "/", gene, "-", condition, "-GWAS-LD.ld"), sep = "") %>% 
    dplyr::rename(
        variant_id = SNP_B,
        r2_lead = R2
    ) %>% 
    select(variant_id, r2_lead)

############
# Prep for plotting
############
# Combine
print("..Putting all together")
both = eqtl %>% 
    select(variant_id, pval_nominal, pos) %>%
    left_join(eqtl_ld) %>% 
    rowwise() %>% 
    mutate(
        type = "eQTL",
        is_lead = variant_id == !!eqtl_lead
    ) %>% 
    bind_rows(
        gwas %>% 
            dplyr::rename(
                variant_id = RSid,
                pval_nominal = pval,
                pos = Pos
            ) %>% 
            select(variant_id, pval_nominal, pos) %>%
            left_join(gwas_ld) %>% 
            mutate(
                type = "GWAS",
                is_lead = variant_id == !!gwas_lead
            )
    ) %>% 
    mutate(logpval = -log10(pval_nominal))

# Get the formatted condition name
condition_tissue = tail(unlist(strsplit(condition, "_")), 1)
if(condition_tissue == "ti"){
    nicetissue = "TI"
} 
if(condition_tissue == "r"){
    nicetissue = "Rectum"
} 
if(condition_tissue == "blood"){
    nicetissue = "Blood"
} 
if(condition_tissue == "ct"){
    nicetissue = "Cross-site"
} 
condition_no_tissue = gsub(paste0("_", condition_tissue), "", condition)
if(grepl("_", condition_no_tissue)){
    nicename = paste0(read.csv(annotation_file) %>% 
                        mutate(label_new = gsub("_", "\\ ", JAMBOREE_ANNOTATION)) %>% 
                        filter(leiden == !!condition_no_tissue) %>% 
                        pull(label_new),
                    " (", nicetissue, ")")
} else {
    nicename = paste0(condition_no_tissue, " major population", " (", nicetissue, ")")
}

# Add new line if too long
wrap_nicename <- function(txt, width = 28) {
  if (nchar(txt) <= width) return(txt)
  
  # find spaces before the cutoff
  spaces <- gregexpr(" ", txt)[[1]]
  cut_pos <- max(spaces[spaces <= width])
  
  # insert line break
  paste0(
    substr(txt, 1, cut_pos - 1), "\n",
    substr(txt, cut_pos + 1, nchar(txt))
  )
}
nicename = wrap_nicename(nicename) # Put on two seperate lines if too long

# Annotate variants in the other clumps by their minimum res
minres = read.delim(eqtl_minres_map) %>% 
  filter(grepl(!!gene, phenotype_clump_index)) %>% 
  rowwise() %>%
  mutate(
    qtl_clump_index = unlist(strsplit(phenotype_clump_index, "\\-"))[c(F,T)]
  ) %>%
  select(qtl_clump_index, annotation_type)

min_var_per_clump = read.delim(clump_file) %>% 
    filter(phenotype_id == !!gene) %>% 
    left_join(minres) %>% 
    left_join(eqtl) %>% 
    group_by(qtl_clump_index) %>% 
    slice_min(pval_nominal, with_ties=F) %>%
    ungroup()

both = both %>% 
  left_join(
    min_var_per_clump %>% 
      select(variant_id, annotation_type)
    )
  
# Make sure annotation_type for GWAS is NA
both$annotation_type = ifelse(both$type == "GWAS", NA, both$annotation_type)
both$annotation_type = factor(both$annotation_type, levels=c("All Cells", "Major population", "Cell type", "ieQTL"))
both$index = !is.na(both$annotation_type)


# Remove GWAS variants not in QTL
both = both %>%
    filter(!is.na(r2_lead))

# Define annotations
coloc_stat = fread(paste0(coloc_dir, gwas_trait, ".gz")) %>%
    filter(
        phenotype_id == !!gene,
        condition_name == paste0("dMean__", condition, "_all")
    ) %>% 
    pull(PP.H4.abf)

qtllab <- paste0(nicename, "\n", symbol)
num_qtllab_lines = length(gregexpr("\n", qtllab)[[1]])
qtllab_just = 0.90-(0.1*num_qtllab_lines)
gwaslab = paste0(gwas_trait, " GWAS\nPP.H4=", sprintf("%.2f", coloc_stat))

# Tidy the temp files
print("..Cleaning temp files")
system(sprintf('rm %s ', paste0(tempout, "/", gene, "*")))

############
# Plot
############
print("..Plotting")
# Define fname and plot
if(plot_minres){
  outf=paste0(out, "/", gene, "-", symbol, "-", condition, "-", gwas_trait, "_manhattan-", "eQTL_minres", ".png")

  # Plot
  scatter_plot = ggplot(both, aes(x = pos, y = logpval)) +
  geom_point(
    aes(fill = r2_lead, shape = is_lead),
    color = "transparent", size = 3, stroke = 1
  ) +
  geom_point(
      data = subset(both, is_lead),
      aes(fill = r2_lead, shape = is_lead),
      color = "black", size = 3, stroke = 1.5
    ) +
  scale_fill_gradient(low = "lightgrey", high = "darkorange2") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none")

  if(plot_minres_lines){
    vline_data <- both %>% 
      filter(!is.na(annotation_type)) %>%
      select(pos, annotation_type) %>%
      distinct()

    scatter_plot = scatter_plot + 
      geom_vline(data = vline_data, 
             aes(xintercept = pos, color = annotation_type), linewidth=1, lty="dashed") +
      scale_color_manual(values=annot.class.palette, name = "Granularity/type")
  } else {
    scatter_plot = scatter_plot + 
      geom_point(
        data = subset(both, index),
        aes(fill = r2_lead, shape = is_lead, color = annotation_type),
        size = 3, stroke = 1.5
      ) +
      scale_color_manual(values=annot.class.palette, name = "Granularity/type") +
      guides(color = guide_legend(override.aes = list(shape = 21, fill = "white", size = 3)))
  }

  scatter_plot = scatter_plot +
    geom_point(
      data = subset(both, is_lead),
      aes(fill = r2_lead, shape = is_lead),
      color = "black", size = 3, stroke = 1.5
    ) # Make sure the lead is black

} else {
  outf=paste0(out, "/", gene, "-", symbol, "-", condition, "-", gwas_trait, "_manhattan.png")

  # Plot
  scatter_plot = ggplot(both, aes(x = pos, y = logpval)) +
    geom_point(
      aes(fill = r2_lead, shape = is_lead),
      color = "transparent", size = 3, stroke = 1
    ) +
    geom_point(
      data = subset(both, index),
      aes(fill = r2_lead, shape = is_lead),
      color = "darkorange", size = 3, stroke = 1
    ) +
    geom_point(
      data = subset(both, is_lead),
      aes(fill = r2_lead, shape = is_lead),
      color = "black", size = 3, stroke = 1
    ) +
    scale_fill_viridis_c() +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") 
}

scatter_plot = scatter_plot + 
  facet_grid(type ~ ., labeller = label_value, scales = "free_y") +
  theme_classic() +
  theme(
    strip.text.y = element_blank(),
    strip.background = element_blank()
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x / 1e6),
    limits = c(min,max)
  ) + 
  geom_text(
    data = data.frame(
      type = unique(both$type),
      pos = min(both$pos),
      logpval = c(
        max(both %>% filter(type == "eQTL") %>% pull(logpval)) * qtllab_just, 
        max(both %>% filter(type == "GWAS") %>% pull(logpval) * 0.85)
        ),
      label = c(qtllab, gwaslab)
    ),
    aes(x = pos, y = logpval, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 0, size=4.5
  ) + 
  labs(
    x=paste0("Chromosome ", chrnum, " position (Mb)"),
    y=expression(-log[10](p-value)),
    fill = "r²"
  ) + 
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=13),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  )

# Define the gene packing function
pack_intervals <- function(df) {
  # Sort genes by start
  df <- df[order(df$start), ]
  rows <- list()
  df$row <- NA_integer_
  
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    for (r in seq_along(rows)) {
      # Check overlap with the last interval on this row
      last_end <- rows[[r]]
      if (df$start[i] > last_end) {
        df$row[i] <- r
        rows[[r]] <- df$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      # Start a new row
      r <- length(rows) + 1
      df$row[i] <- r
      rows[[r]] <- df$end[i]
    }
  }
  df
}

# get gene pos
gene_df <- read.delim(gene_pos_f) %>%
  mutate(
    tss   = ifelse(strand == "+", start, end),
    other = ifelse(strand == "+", end, start)
  ) %>% 
  filter(
    chromosome == chrnum,
    tss > min,
    tss < max,
    end > min,
    end < max
  ) %>%
  pack_intervals() %>% 
  rowwise() %>%
  mutate(
    want = feature_id == !!gene,
    label = ifelse(want, !!symbol, "")
  )

# Plot the gene positions
gene_plot <- ggplot(gene_df) +
  geom_segment(
    aes(x = tss, xend = other,
        y = row, yend = row,
        color = want),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    size = 1.2
  ) +
  geom_text(data = gene_df,
            aes(x = (start + end)/2, y = row+0.4, label = label),
            inherit.aes = FALSE, hjust = 0, size = 4, color="darkorange2") +
  scale_color_manual(
    values = c(`TRUE` = "darkorange2", `FALSE` = "lightgrey"),
    guide = "none"
  ) + 
  labs(x = paste0("Chromosome ", chrnum, " position (Mb)"),
       y = NULL) +
  ylim(c(0.5,max(gene_df$row)+0.5)) + 
  theme_classic() +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x / 1e6),
    limits = c(min,max)
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=12),
    axis.title.x = element_text(size=13)
  )

if(plot_minres_lines){
  gene_plot <- gene_plot +
    geom_vline(data = vline_data,
               aes(xintercept = pos, color = annotation_type),
               linewidth = 1, lty = "dashed") +
    scale_color_manual(values = annot.class.palette, guide = "none") +
    scale_color_manual(
      values = c(`TRUE` = "darkorange2", `FALSE` = "black"),
      guide = "none"
    ) # Hide legend on gene_plot
}

# Construct the NeQTL per res fig
plot_variants = both %>% 
  filter(!is.na(annotation_type)) %>%
  distinct(variant_id, pos, annotation_type) %>% 
  left_join(
    read.delim(clump_file) %>%
      filter(phenotype_id == !!gene)
  ) 

eqtl_count = read.delim(eQTL_count_map) %>%
      filter(grepl(!!gene, phenotype_clump_index)) %>%
      rename(qtl_clump_index = "phenotype_clump_index") %>%
      mutate(qtl_clump_index = unlist(strsplit(qtl_clump_index, "\\-"))[c(F,T)]) %>%
      left_join(
        plot_variants %>% 
          select(qtl_clump_index, pos)
      ) %>%
      mutate(annotation_type = factor(annotation_type, levels = c("Cell type", "Major population", "All Cells", "ieQTL")))
  
eqtl_count_plot = ggplot(eqtl_count, aes(x = pos, y = annotation_type)) + 
  geom_tile(aes(fill = count), width = 15000, height = 1, color="black") +
  scale_fill_viridis_c(limits=c(0, max(eqtl_count$count)), name="Number of\nannotations") +  
  theme_classic() +
  labs(x = "Position", y = "", fill = "Count") + 
  scale_x_continuous(limits = c(min,max)) + 
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12),
    axis.title.y = element_blank()
  )

# Combine
if(plot_eqtl_count){
  eqtl_count_plot / scatter_plot / gene_plot +
  plot_layout(heights = c(0.8, 3, 0.5), guides = "collect")   
} else {
  scatter_plot / gene_plot +
  plot_layout(heights = c(3, 0.5))   
}


ggsave(outf, width = 10, height = 7.5)
print("..DONE!")

