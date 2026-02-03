# IBDVerse-sc-eQTL-code

Contained in this repository are the code and instructions to reproduce the single-cell eQTL results present in the preprint: Cell-type-resolved genetic regulatory variation shapes inflammatory bowel disease risk (Alegbe, Harris et al., 2025).

## Performing quality control and normalization of single-cell RNA-sequencing data

We used the yascp pipeline to perform initial data processing, quality control and normalization of our single-cell RNA-sequencing data. The yascp package is available at: https://github.com/wtsi-hgi/yascp. We used version 1.3.0 of yascp for our analyses.

The config files can be found within this repository in the `configs/` folder with names starting with `yascp`. The config files specify the input data, quality control thresholds, and parameters used for normalization.

To run yascp, please follow the instructions in the yascp repository.

## Quality control and clustering of single-cell RNA-sequencing data

Code used to perform all quality control and clustering of scRNAseq data is available at https://github.com/andersonlab/atlassing. 

## Mapping single-cell eQTLs

We used the QTLight pipeline to pseudobulk our single-cell RNA-sequencing data and map eQTLs with TensorQTL. The package QTLight is available at: https://github.com/wtsi-hgi/QTLight. We used version 1.4.0 of QTLight for our analyses.

The config files can be found within this repository in the `configs/` folder with names starting with `qtlmap`. The config files specify the input data, covariates, and parameters used for mapping eQTLs. 

To run QTLight, please follow the instructions in the QTLight repository.

## Post-processing of eQTL results

The code for testing eQTL replication rates is available at https://github.com/andersonlab/pi1_pairs, for colocalisation at https://github.com/andersonlab/snakemake_colocalisation_sceqtl, and for LD clumping at https://github.com/andersonlab/ld_clump_tqtl.

The remaining code for generating plots and performing additional analyses is available in this repository. If anything is unclear or missing, please open an issue or contact Tobi Alegbe (Tobi1Kenobi) or Bradley Harris (BradleyH017) directly.



