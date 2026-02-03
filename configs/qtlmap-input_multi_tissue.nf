params{
    method= 'single_cell' //or a [bulk | single_cell] (if single cell used the *phenotype_file* is a h5ad file)
    input_vcf ='batch5_imputed.vcf.gz' // Imputed genotypes
    genotype_phenotype_mapping_file = '' // annotation file containing Genotype | RNA and | Condiion
    annotation_file = './IBDVerse-sc-eQTL-code/data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt' //assets file that has start and end positions for the genes, this one is using hg38
    phenotype_file = 'cross_tissue-analysis_anndata.h5ad' //this should point to h5ad file in a single cell experiments.
    aggregation_columns='unannotated,predicted_category,predicted_labels,unannotated_tissue,predicted_category_tissue,predicted_labels_tissue' // For base
    // aggregation_columns='unannotated_tissue,predicted_category_tissue,predicted_labels_tissue' // For interaction
    gt_id_column='Genotyping_ID' //for the scrna h5ad defines which column has individual level id (corresponding to the VCF file)
    sample_column='Genotyping_ID' // for the scrna h5ad defines which column has the sample name (can be identical to gt_id_column if sample=individual)
    sample_covariates='' //covariates to be included in the model - LEAVE BLANK, DOES NOTHING AT PRESENT
}