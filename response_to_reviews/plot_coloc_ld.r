# Bradley August 2025
# Plot coloc plots, colour points by LD
# module load HGI/softpack/users/bh18/Colocalisation_analysis/2
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000184384" "MAML2" "Myeloid_4_ti" "UC"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000172575" "RASGRP1" "Colonocyte_ct" "CD" 
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000136997" "MYC" "Epithelial_12_ct" "CD"  
#
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000108175" "ZMIZ1" "T_r" "CD" 
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000162613" "FUBP1" "Colonocyte_r" "CD" 
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000101311" "FERMT1" "Secretory_ct" "IBD"  
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000137806" "NDUFAF1" "Secretory_ti" "UC" 
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000143801" "PSEN2" "Epithelial_3_r" "UC"  
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000143801" "PSEN2" "Epithelial_1_ti" "UC" #  (out of interest)

# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000013561" "RNF14" "Epithelial_1_ti" "CD" "ses_cd_binned"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000115232" "ITGA4" "Myeloid_blood" "IBD"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000096968" "JAK2" "Epithelial_23_ct" "IBD"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000163513" "TGFBR2" "Epithelial_5_ct" "IBD"

# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000163513" "MOSPD3" "T_6_ct" "UC"
# Rscript response_to_reviews/plot_coloc_ld.r "ENSG00000163513" "MOSPD3" "Enterocyte_ct" "UC"

################
# Preamble
################
# Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(scales))
# Truncate magma at orange
pal <- viridisLite::magma(256)      # full palette
pal_trunc <- pal[1:200]   


# Options
args = commandArgs(trailingOnly=TRUE)
gene = args[1] # gene="ENSG00000162613"
symbol = args[2] # symbol = "FUBP1"
condition = args[3] # condition="Colonocyte_r"
gwas_trait = args[4] # gwas_trait="CD"
if(length(args) == 5){
  interaction = args[5] # interaction = "ses_cd_binned"
} else {
  interaction = NULL
}
print(paste0("~~~~~~ Plotting coloc between ", symbol, " and ", gwas_trait, " in ", condition))

# Hard code paths and options
sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
sumstats.interaction.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_12-multi_tissue_interaction_results/TensorQTL_eQTLS/'
gwas_dir <- "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/gwas_sumstats_final/"
gwas_master_list = paste0(gwas_dir, "gwas_files.txt")
clump_file = "data/clumped_all.txt.gz"
window = 5e5 # Centred on eQTL lead, in either direction
plink_prefix = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/plink/imputed_chr" # Location of plink genotyping files
annotation_file = "data/all_IBDverse_annotation_mastersheet.csv"
out = "response_to_reviews/coloc_plots"
gene_pos_f = "data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
coloc_dir = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/'
int_coloc_dir = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed/'
tempout = "temp/coloc"
eqtl_minres_map = "eqtl_out/eQTL_minres_map.txt.gz"
eQTL_count_map = "eqtl_out/eQTL_count_map.txt.gz"
for(dir in c(out, tempout)){
    if(!file.exists(dir)){
        dir.create(dir)
    }
}
annot.class.palette <- c('All Cells'='#edf8b1',
                         'Major pop.'='#7fcdbb',
                         'Cell type'='#2c7fb8')

interaction_name_convert = data.frame(
  raw = c("age_binned", "ses_cd_binned", "ses_inflamed", "sex", "smoking_ever"),
  formatted = c("Age", "Inflam. status\n(binned)", "Inflam. status", "Sex", "Smoker")
)

# Define variables based on whether interaction or base
if(is.null(interaction)){
  eqtl_sumstatpath = paste0(sumstats.all.basedir, "dMean__", condition, "_all/OPTIM_pcs/base_output/base/")
  eqtlf = paste0(eqtl_sumstatpath, "Cis_eqtls_qval.tsv")
  nom_pval_col = "pval_nominal"
  coloc_path = paste0(coloc_dir, gwas_trait, ".gz")
  condition_name_format = paste0("dMean__", condition, "_all")
  outf=paste0(out, "/", gene, "-", symbol, "-", condition, "-", gwas_trait, "_manhattan-", "eQTL_minres", ".png")
} else {
  eqtl_sumstatpath = paste0(sumstats.interaction.basedir, "dMean__", condition, "_all/OPTIM_pcs/interaction_output/", interaction, "/")
  eqtlf = paste0(eqtl_sumstatpath, "cis_inter1.cis_qtl_top_assoc.txt.gz")
  nom_pval_col = "pval_gi"
  coloc_path = paste0(int_coloc_dir, gwas_trait, ".gz")
  condition_name_format = paste0("dMean__", condition, "_all__", interaction)
  niceintname = interaction_name_convert %>% 
    filter(raw == !!interaction) %>% 
    pull(formatted)
  annot.class.palette = c(annot.class.palette, 'ieQTL' = '#f8b1f1ff')
  outf=paste0(out, "/", gene, "-", symbol, "-", condition, "-", interaction, "-", gwas_trait, "_manhattan-", "eQTL_minres", ".png")
}


############
# Get gwas and QTL sumstats
############
## eQTL
# Lead
print("..Loading the eQTL sumstats")
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

if(is.null(interaction)){
  eqtlnom = paste0(eqtl_sumstatpath, "cis_nominal1.cis_qtl_pairs.", chrnum, ".tsv")
} else {
  eqtlnom = paste0(eqtl_sumstatpath, "cis_inter1.cis_qtl_pairs.", chrnum, ".tsv")
}


# Nominal
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
    filter(RSid %in% eqtl$variant_id) %>% 
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
    select(variant_id, !!nom_pval_col, pos) %>%
    rename(pval_nominal = !!nom_pval_col) %>% 
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
    nicename = paste0(condition_no_tissue, " Major Pop.", " (", nicetissue, ")")
}

# Add new line if too long
wrap_nicename <- function(txt, width = 40) {
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
nicename = gsub("CAECAM", "CEACAM", nicename)

# Annotate variants in the other clumps by their minimum res
minres = read.delim(eqtl_minres_map) %>% 
  filter(grepl(!!gene, phenotype_clump_index)) %>% 
  rowwise() %>%
  mutate(
    qtl_clump_index = unlist(strsplit(phenotype_clump_index, "\\-"))[c(F,T)]
  ) %>%
  select(qtl_clump_index, annotation_type)

# Exclude ieQTLs if not working with interaction
if(is.null(interaction)){
  minres = minres %>%
    filter(annotation_type != "ieQTL")
}

min_var_per_clump = read.delim(clump_file) %>% 
    filter(phenotype_id == !!gene) %>% 
    left_join(minres) %>% 
    left_join(eqtl %>% 
      rename(pval_nominal = !!nom_pval_col)) %>% 
    filter(!is.na(pval_nominal)) %>%
    group_by(qtl_clump_index) %>% 
    slice_min(pval_nominal, with_ties=F) %>%
    ungroup()

both = both %>% 
  left_join(
    min_var_per_clump %>% 
      select(variant_id, annotation_type)
    )

# Replace the resolutions to use abbreviation
both$annotation_type = gsub("Major population", "Major pop.", both$annotation_type)

# Make sure annotation_type for GWAS is NA
both$annotation_type = ifelse(both$type == "GWAS", NA, both$annotation_type)
both$annotation_type = factor(both$annotation_type, levels=names(annot.class.palette))
both$index = !is.na(both$annotation_type)


# Remove GWAS variants not in QTL
both = both %>%
    filter(!is.na(r2_lead))

# Define annotations
coloc_stat = fread(coloc_path) %>%
    filter(
        phenotype_id == !!gene,
        condition_name == !!condition_name_format
    ) %>% 
    pull(PP.H4.abf)

if(is.null(interaction)){
  qtllab <- paste0("italic('", symbol, "')~'eQTL'")
} else {
  qtllab <- paste0("italic('", symbol, "')~'ieQTL'")
}
gwaslab = paste0(gwas_trait, " GWAS\nPP.H4=", sprintf("%.2f", coloc_stat))

# Tidy the temp files
print("..Cleaning temp files")
system(sprintf('rm %s ', paste0(tempout, "/", gene, "*")))

# Create flag for the eQTL lead
both = both %>%
  mutate(is_eqtl_lead = type == "eQTL" & is_lead)

# For the plotting of lead with seperate colour, need to discretise these
both = both %>% mutate(
    # Create bins for r2_lead (adjust breaks as needed)
    r2_binned = cut(r2_lead, breaks = seq(0, 1, 0.2), include.lowest = TRUE),
    
    # Create combined variable for coloring
    color_var = ifelse(is_eqtl_lead, 
                      as.character(annotation_type),
                      as.character(r2_binned))
  )

r2_breaks <- seq(0, 1, 0.2)
r2_colors <- pal_trunc[seq(1, length(pal_trunc), length.out = length(r2_breaks)-1)]
names(r2_colors) <- paste0(levels(cut(r2_breaks, breaks = r2_breaks)))
combined_palette <- c(r2_colors, annot.class.palette)

# Format these (make nice and relevel)
names(combined_palette) = gsub("\\(", "", names(combined_palette))
names(combined_palette) = gsub("\\[", "", names(combined_palette))
names(combined_palette) = gsub("\\]", "", names(combined_palette))
names(combined_palette) = gsub("\\,", " - ", names(combined_palette))
both$color_var = gsub("\\(", "", both$color_var)
both$color_var = gsub("\\[", "", both$color_var)
both$color_var = gsub("\\]", "", both$color_var)
both$color_var = gsub("\\,", " - ", both$color_var)
both$color_var = factor(both$color_var, levels=c(c("0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2", names(annot.class.palette))))

# Make sure all anotation type values are plotted:
both$annotation_type <- factor(both$annotation_type, levels = names(annot.class.palette))

# Make a dummy for the colour legend (plot behind the left-most point)
mineqtlpos = both %>%
  filter(type == "eQTL") %>% 
  arrange(pos) %>% 
  head(n=1)

legend_data_col <- data.frame(
  annotation_type = names(annot.class.palette),
  type = "eQTL",
  logpval = mineqtlpos$logpval,
  pos=mineqtlpos$pos
)


# Do the same for the LD
mingwaspos = both %>%
  filter(type == "GWAS") %>% 
  arrange(pos) %>% 
  head(n=1)

legend_data_fill <- data.frame(
  color_var = c("0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"),
  type="GWAS",
  logpval=mingwaspos$logpval,
  pos=mingwaspos$pos
)
legend_data_fill$color_var = factor(legend_data_fill$color_var, levels = c("0.8 - 1", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0 - 0.2"))
# Reorder palette
ldcols = rev(combined_palette[grep("\\-", names(combined_palette))])
rescols = combined_palette[-grep("\\-", names(combined_palette))]
combined_palette = c(ldcols, rescols)

############
# Plot
############
print("..Plotting")

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
    label = ifelse(want, !!symbol, ""),
    label = ifelse(label != "", paste0("italic('", label, "')"), "")
  )

# Plot the gene positions
gene_plot <- ggplot(gene_df) +
  geom_segment(
    aes(x = tss, xend = other,
        y = row, yend = row,
        color = want),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    size = 1.7
  ) +
  geom_text(data = gene_df,
            aes(x = (start + end)/2, y = row+0.4, label = label),
            inherit.aes = FALSE, hjust = 0, size = 6, color="darkorange2", parse=T) +
  scale_color_manual(
    values = c(`TRUE` = "darkorange2", `FALSE` = "black"),
    guide = "none"
  ) + 
  labs(x = paste0("Chromosome ", chrnum, " position (Mb)"),
       y = NULL) +
  ylim(c(0.5,max(gene_df$row)+0.5)) + 
  theme_classic() +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x / 1e6),
    limits = c(min,max),
    expand = c(0.02, 0.02)
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=18),
    axis.title.x = element_text(size=20)
  )

# Make top plot (eQTL)
topdat = both %>% filter(type == "eQTL")
top = ggplot(topdat, aes(x = pos, y = logpval)) +
  geom_point(data = legend_data_col, # Plot the dummy first
          aes(pos, logpval, colour = annotation_type)) +
  geom_point( # Plot the points for non-leads, coloured by LD
    data = topdat %>% filter(!is_eqtl_lead),
    aes(fill = color_var, shape = is_lead),
    color = "transparent", size = 3, stroke = 1
  ) +
  geom_point( # Plot the leads with black outline
    data = topdat %>% filter(is_lead),
    aes(fill = color_var, shape = is_lead),
    color = "black", size = 4, stroke = 1.5
  ) +
  scale_fill_manual(values = combined_palette,
    breaks = names(combined_palette)[grep("\\-", names(combined_palette))],
  ) + 
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") + 
  geom_point(
      data = topdat %>% filter(index, !is_eqtl_lead),
      aes(fill = color_var, shape = is_lead, color = annotation_type),
      size = 3, stroke = 2
    ) +
  scale_color_manual(
    values=annot.class.palette, 
    name = "Resolution", 
    breaks = names(annot.class.palette),
    drop=F
  ) +
  guides(
    color = "none",
    fill = guide_legend(
      title = "r²",
      override.aes = list(
        shape = 21, 
        color = "transparent",
        size = 5
      ),
      title.theme = element_text(size = 19),
      label.theme = element_text(size = 16)
    )
  ) + 
  theme_classic() + 
  theme(
      legend.position = c(0.95, 1.025),     # top-right inside panel
      legend.justification = c("right","top"),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size=19),
      axis.title.y = element_text(size=21),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
      axis.line.x = element_blank(),
      title = element_text(size=18)
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x / 1e6),
    limits = c(min,max),
    expand = c(0.02, 0.02)
  ) + 
  geom_text(
    data = data.frame(
      type = "eQTL",
      pos = both %>% filter(type == "eQTL") %>% pull(pos) %>% min(),
      logpval = both %>% filter(type == "eQTL") %>% pull(logpval) %>% max()*0.95,
      label = c(qtllab)
    ),
    aes(x = pos, y = logpval, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 0, size=7, parse=TRUE
  ) + 
  labs(
    x=paste0("Chromosome ", chrnum, " position (Mb)"),
    y=expression(-log[10](p-value)),
    fill = "r²",
    title = nicename
  )

# Plot extra annotation 
if(!is.null(interaction)){
  top = top +  
    geom_text(
      data = data.frame(
        type = "eQTL",
        pos = both %>% filter(type == "eQTL") %>% pull(pos) %>% min(),
        logpval = both %>% filter(type == "eQTL") %>% pull(logpval) %>% max()*0.68,
        label = c(niceintname)
      ),
      aes(x = pos, y = logpval, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 0, size=7
    ) 
}

midat = both %>% filter(type == "GWAS")
middle = ggplot(midat, aes(x = pos, y = logpval)) +
  geom_point(data = legend_data_fill, # Plot the dummy first
          aes(pos, logpval, fill = color_var)) +
  geom_point( # Plot the points for non-leads, coloured by LD
    data = midat %>% filter(!is_eqtl_lead),
    aes(fill = color_var, shape = is_lead),
    color = "transparent", size = 3, stroke = 1
  ) +
  geom_point( # Plot the leads with black outline
    data = midat %>% filter(is_lead),
    aes(fill = color_var, shape = is_lead),
    color = "black", size = 4, stroke = 1.5
  ) +
  scale_fill_manual(values = combined_palette,
  breaks = names(combined_palette)[grep("\\-", names(combined_palette))],
  ) + 
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") + 
  geom_point(
      data = midat %>% filter(index, !is_eqtl_lead),
      aes(fill = color_var, shape = is_lead, color = annotation_type),
      size = 3, stroke = 2
    ) +
  scale_color_manual(
    values=annot.class.palette, 
    name = "Resolution", 
    breaks = names(annot.class.palette),
    drop=F
  ) +
  guides(
    fill = "none",
    color = guide_legend(
      title = "Resolution",
      override.aes = list(
        shape = 21, 
        fill = "white", 
        size = 5,         
        stroke = 2
      ),
      title.theme = element_text(size = 19),      
      label.theme = element_text(size = 16)       
    )
  ) + 
  theme_classic() + 
  theme(
      legend.position = c(0.95, 1.025),     # top-right inside panel
      legend.justification = c("right","top"),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(size=19),
      axis.title.y = element_text(size=21),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  ) +
  scale_x_continuous( 
    labels = function(x) sprintf("%.1f", x / 1e6),
    limits = c(min,max),
    expand = c(0.02, 0.02)
  ) + 
  geom_text(
    data = data.frame(
      type = "GWAS",
      pos = both %>% filter(type == "GWAS") %>% pull(pos) %>% min(),
      logpval = both %>% filter(type == "GWAS") %>% pull(logpval) %>% max() * 0.82,
      label = c(gwaslab)
    ),
    aes(x = pos, y = logpval, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 0, size=7
  ) + 
  labs(
    x=paste0("Chromosome ", chrnum, " position (Mb)"),
    y=expression(-log[10](p-value)),
    fill = "r²"
  )

top / middle / gene_plot +
  plot_layout(axes = "collect", heights = c(1.5, 1.5, 0.5))   

ggsave(outf, width = 8.5, height = 8.5)
print("..DONE!")