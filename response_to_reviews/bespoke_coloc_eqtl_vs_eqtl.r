# Bradley September 2025
# module load HGI/softpack/groups/hgi/ubeR


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(patchwork))


# Define run options
phenotype = "ENSG00000162613"
symbol = "FUBP1"
coloc_clump = "chr1:77984833:C:A" # index variant from the colocalising clump
conditions = c("Myeloid_8_r", "Mesenchymal_2_r", "T_2_blood") # conditions we find the effect in (cluster x tissue)

# Hard code paths for other bits
sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
gene_pos_f = "data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
gene_pos = read.delim(gene_pos_f)
annotation_mastersheet = "data/all_IBDverse_annotation_mastersheet.csv"
clump_file = paste0("data/clumped_all.txt.gz")


outf = paste0("eqtl_out/", symbol, "-", paste(conditions, collapse="."), ".png")

annots = read.csv(annotation_mastersheet)
chrnum = gene_pos %>%
    filter(feature_id == !!phenotype) %>% 
    pull(chromosome) %>%
    head(1)


# Load the nominal eQTL effects for this gene in each condition
nominal = lapply(conditions, function(condition){
    path = paste0(sumstats.all.basedir, "dMean__", condition, "_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.", chrnum, ".tsv")
    return(fread(path) %>% 
        filter(phenotype_id == !!phenotype) %>% 
        mutate(name = !!condition))
})
names(nominal) = conditions

# For each pairs of conditions, perform coloc
comparisons = combn(as.character(conditions), 2, simplify = FALSE)
allres = NULL
for(pair in comparisons){
    # Set up dataset
    D1=list(
        beta = nominal[[pair[1]]]$slope,
        varbeta = nominal[[pair[1]]]$slope_se^2, # varbeta is SE squared
        snp = nominal[[pair[1]]]$variant_id,
        position = unlist(strsplit(nominal[[pair[1]]]$variant_id, "\\:"))[c(F,T,F,F)],
        sdY = 1, 
        type = "quant"
    )

    D2=list(
        beta = nominal[[pair[2]]]$slope,
        varbeta = nominal[[pair[2]]]$slope_se^2,
        snp = nominal[[pair[2]]]$variant_id,
        position = unlist(strsplit(nominal[[pair[2]]]$variant_id, "\\:"))[c(F,T,F,F)],
        sdY = 1,
        type = "quant"
    )

    res = coloc.abf(dataset1=D1, dataset2=D2)$summary %>% t() %>% as.data.frame() # Default prior
    res$cond1 = pair[1]
    res$cond2 = pair[2]
    rownames(res) = NULL

    if(is.null(allres)){
        allres = res
    } else {
        allres = rbind(allres, res)
    }
}


# Plot the associations at these loci + for all other leads - show these are truly not the same signal
clumps = clumps = read.delim(clump_file) %>% 
    filter(phenotype_id == !!phenotype) %>%
    mutate(pos = as.numeric(unlist(strsplit(variant_id, "\\:"))[c(F,T,F,F)]))


plot_dat = do.call(rbind, nominal) %>% 
    mutate(
        pos = as.numeric(unlist(strsplit(variant_id, "\\:"))[c(F,T,F,F)]),
        int_var = variant_id == variant_of_interest,
        leiden = gsub("_ti", "", name),
        leiden = gsub("_r", "", leiden),
        leiden = gsub("_blood", "", leiden),
    ) %>% 
    left_join(
        annots %>%
            select(leiden, JAMBOREE_ANNOTATION) %>%
            mutate(
                JAMBOREE_ANNOTATION = gsub("\\_", " ", JAMBOREE_ANNOTATION)) 
    ) %>%
    mutate(
        tissue = ifelse(grepl("_r", name), "Rectum", ""),
        tissue = ifelse(grepl("_ti", name), "TI", tissue),
        tissue = ifelse(grepl("_blood", name), "Blood", tissue),
        nicename = paste0(JAMBOREE_ANNOTATION, "-", tissue)
    ) %>%
    left_join(
        clumps %>% 
            select(variant_id, qtl_clump_index)
    ) %>%
    mutate(
        is_coloc_clump = qtl_clump_index == !!coloc_clump
    )

maxpos = min(max(clumps$pos+5e5), max(plot_dat$pos))
minpos = max(min(clumps$pos-5e5), min(plot_dat$pos))

plot_dat = plot_dat %>% 
    filter(
        pos > minpos,
        pos < maxpos
    )

# Get min pvalue per clump, annotate these
min_per_clump = plot_dat %>% 
    group_by(nicename, qtl_clump_index) %>% 
    filter(!is.na(qtl_clump_index)) %>%
    slice_min(pval_nominal, with_ties=F) %>%
    select(nicename, variant_id, qtl_clump_index) %>%
    mutate(is_min = T)

plot_dat = plot_dat %>% 
    left_join(min_per_clump) %>%
    mutate(
        assoc_type = ifelse(variant_id == coloc_clump, "coloc", ""),
        assoc_type = ifelse(is_min == T & is_coloc_clump == F, "other_qtl", assoc_type),
        assoc_type = ifelse(is.na(assoc_type), "none", assoc_type)
    )

ggplot(plot_dat, aes(x=pos, y=-log10(pval_nominal))) + 
    geom_point(
        aes(fill = assoc_type, shape = assoc_type),
        color = "transparent", size = 3, stroke = 1
    ) + 
    geom_point(
        data = subset(plot_dat, assoc_type == "other_qtl"),
        aes(fill = assoc_type, shape = assoc_type),
        color = "transparent", size = 3, stroke = 1.5
    ) +
    geom_point(
        data = subset(plot_dat, assoc_type == "coloc"),
        aes(fill = assoc_type, shape = assoc_type),
        color = "black", size = 3, stroke = 1.5
    ) +
    scale_shape_manual(values = c(none = 21, other_qtl = 21, coloc = 23), guide = "none") + 
    scale_fill_manual(values = c(none = "lightgrey", other_qtl = "black", coloc = "darkorange2"), guide = "none") + 
    facet_grid(name ~ ., labeller = label_value, scales = "free_y") +
    theme_classic() +
    theme(
        strip.text.y = element_blank(),
        strip.background = element_blank()
    ) +
    scale_x_continuous(
        labels = function(x) sprintf("%.1f", x / 1e6)
    ) + 
    geom_text(
    data = data.frame(
        name = unique(plot_dat$name),  # Add the faceting variable
        type = unique(plot_dat$nicename),
        pos = maxpos,
        logpval = c(
            plot_dat %>% 
                group_by(nicename) %>%
                summarise(maxp = max(-log10(pval_nominal))*0.85) %>%
                arrange(nicename) %>%
                pull(maxp)
        )
    ),
    aes(x = pos, y = logpval, label = type),
    inherit.aes = FALSE, hjust = 1, vjust = 0, size=4.5
    ) + 
    labs(
        title = paste0(symbol, ". Max PP.H4=", signif(max(allres$PP.H4.abf), 2)),
        x=paste0("Chromosome ", chrnum, " position (Mb)"),
        y=expression(-log[10](p-value))
    )

ggsave(outf, width = 8, height = 7)
print("..DONE!")
