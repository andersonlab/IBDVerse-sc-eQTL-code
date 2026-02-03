# Bradley July 2025
# module load $scvi

# Packages
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.pyplot import rc_context
from anndata.experimental import read_elem
from h5py import File

########
# Format outdir
########
outdir="eqtl_out"
sc.settings.figdir=outdir
sc.settings.verbosity = 3    
sc.logging.print_header()
sc.settings.set_figure_params(dpi=500, facecolor='white', format="png")


#################
# Checking sequencing depth of healthy vs control samples
#################
# Load data (obs only)
h5ad = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/adata_PCAd_batched_umap_add_expression.h5ad"
f2 = File(h5ad, 'r')
obs = read_elem(f2['obs'])
# Subset to ti
obs = obs[obs['tissue'] == "ti"] # 913569 cells

# Append the number of cells per sample
obs['ncells'] = obs['sanger_sample_id'].map(obs['sanger_sample_id'].value_counts())

# Extract what we want to save
df = obs[["sanger_sample_id", "tissue", "disease_status", "Median_nCounts", "ncells"]].reset_index(drop=True).drop_duplicates()
df.to_csv(f"{outdir}/per_sample_data_ti_atlas.csv")


##################
# Making supp metadata file
##################
h5ad="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/atlassing/results/objects/from_irods/celltypist_prediction.h5ad"
f2 = File(h5ad, 'r')
obs = read_elem(f2['obs']) # >4M cells

# Get the combination of sample, tissue and Genotyping ID
meta = obs[['sanger_sample_id','Genotyping_ID','tissue', 'disease_status', 'inflammation_status', 'sex' 'age']].reset_index(drop=True).drop_duplicates()

# Match this with the meta for other variables 
other = pd.read_csv("/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/interaction_files/2024_12_27-multi_tissue/eQTL_interactions.tsv", sep ="\t")
meta = meta.merge(other, how="left")