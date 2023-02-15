
########################################################################
## Cell Deconvolution with Nowakowski_2017 - Prepare signature matrix ##
########################################################################

# Samantha O'Connor
# School of Biological Health Systems and Engineering
# Arizona State University
# Plaisier Lab
# saoconn1@asu.edu
# Last updated: September 14, 2022

#############################################

# docker pull cplaisier/scrna_seq_velocity_6_7_2021
# docker run -it -v '/home/soconnor/old_home/scGlioma/Nowakowski/:/files' cplaisier/scrna_seq_velocity_6_7_2021

#------------------------------------------------------------------
# Set up section / load packages
#------------------------------------------------------------------

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom
from scipy.stats import spearmanr
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Set directories
tag = 'Nowakowski'
dir = 'final_output'

#------------------------------------------------------------------
# Prepare signature matrix (Nowakowski combined cell type medioids)
#------------------------------------------------------------------

# Import marker genes for each cell type
mgenes = pd.read_csv(dir+'/'+tag+'_scTransform_Markers_together_celltypes_combined.csv')
df = mgenes.loc[mgenes['p_val_adj'] <= 0.05]
mgenes = df

# Import data (made from prepare_ref_cellDeconvolution.R)
data = sc.read_loom(dir+'/'+tag+'_normalized_celltypes_combined.loom')
data.shape #(3047, 31506)
data.obs['cell_type_combined'].value_counts()

# Filter var_names to marker genes
data_subset = data[:,data.var_names[data.var_names.isin(mgenes['gene'])]]
data_subset.shape #(3047, 9191)
sc.pp.filter_cells(data_subset, min_genes=1)
genes1 = data_subset.shape[1]

# Generate medioids (argmin)
medioids = pd.DataFrame(index=data_subset.var_names, columns=list(data_subset.obs['cell_type_combined'].unique()))
for clust in list(data_subset.obs['cell_type_combined'].unique()):
    print(clust)
    cells = data_subset[data_subset.obs['cell_type_combined']==clust]
    if len(cells.obs) > 1:
        ind = np.argmin(cells.X.sum(axis=0))
        medioids[clust] = cells[ind].X.todense().T

# Save out medioids to csv / tsv with specific format for Cibertsortx import
medioids.to_csv(dir+'/'+tag+'_mgenes_medioids.csv')
medioids = pd.read_csv(dir+'/'+tag+'_mgenes_medioids.csv', index_col=0)
medioids.index.name = 'GeneSymbols'
medioids.to_csv(dir+'/'+tag+'_mgenes_medioids.tsv', sep="\t") # file to compare data to in cibertsortx

# Drop glycolysis cell type and sort for formatting
medioids.drop('Glyc', inplace=True, axis=1)
medioids = medioids[['Neuronal_Progenitor', 'Excitatory', 'Inhibitory', 'Glial_Progenitor', 'Astrocyte', 'Microglia', 'Endothelial', 'Mural', 'Choroid', 'U1']]
medioids.to_csv(dir+'/'+tag+'_mgenes_medioids_noGlycolysis.tsv', sep="\t")


#------------------------------------------------------------------
# Prepare data to deconvolute (mixture file)
#------------------------------------------------------------------

# Import normalized RNA-seq data
data = pd.read_csv("ALS_Control_DESeq2Norm_CIBERSORT_5ksymbols.csv", index_col = 0) #4912
data.index.name = 'GeneSymbol'
data.to_csv(dir+"/ALS_Control_DESeq2_MedianofRatios_Counts_5ksymbols.tsv", sep="\t")
