
###########################################################################
## Cell Deconvolution with Nowakowski_2017 - Plot ALS deconvolution plot ##
###########################################################################

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
fold = 'cibersortx'


#------------------------------------------------------------------
# Prepare data for plotting
#------------------------------------------------------------------

# Load ALS phenotype data
pheno = pd.read_csv("AllSubjects_Phenotype.csv", index_col = 0)

# Load cibersortx data and prepare data for plotting
decon = pd.read_csv(fold+"/CIBERSORTx_091322_Results.csv")
decon2 = decon.merge(pheno[['Subtype']], left_on = 'Mixture', right_on = 'Subject').set_index('Subtype')
decon2.to_csv(dir+"/CIBERSORTx_ALS5000_with_subtype.csv")
decon2.drop('Mixture', axis=1, inplace=True)
deconDF = decon2.T
celltypes1 = ['Glial_Progenitor','Excitatory','Neuronal_Progenitor','Inhibitory','Mural','Endothelial','Astrocyte', 'Microglia', 'Glyc', 'U1', 'Choroid']
deconDF = deconDF.loc[[True if i in celltypes1 else False for i in list(deconDF.index)]]
deconDF2 = deconDF.T
# Only specific subtypes
#subtypes1 = ['OX', 'TE', 'GLIA']
#deconDF2 = deconDF2.loc[[True if i in subtypes1 else False for i in list(deconDF2.index)]]
deconDF2.reset_index(inplace=True)
deconDF2 = deconDF2[['Subtype', 'Neuronal_Progenitor', 'Excitatory', 'Inhibitory', 'Glial_Progenitor', 'Astrocyte', 'Microglia', 'Endothelial', 'Mural', 'Choroid', 'U1']]
df_long = pd.melt(deconDF2, id_vars=['Subtype'])
s1 = df_long.Subtype.unique().tolist()
# Only specific subtypes
#s1_colors=["#b22222","#000080","#daa520"]


#------------------------------------------------------------------
# Generate barplot
#------------------------------------------------------------------

s1_colors=["#b22222","#000080","#ffc125","#FF5F1F", "#3db270"]
color_dict = dict(zip(s1, s1_colors))
plt.figure(figsize=(15,10))
sns.set_theme(style="whitegrid")
sns.boxplot(x='variable', y='value', hue = 'Subtype', linewidth=0.7, fliersize=2, palette = color_dict, data=df_long)
sns.stripplot(x='variable', y='value', hue = 'Subtype', dodge=True, size=2, palette = color_dict, data=df_long)
plt.xticks(rotation = 30)
plt.ylim([0.00, 0.50])
# Only specific subtypes
#plt.savefig(dir+'/'+tag+'_boxplot_swarmplot_5000_genes_specific_subtypes.pdf')
plt.savefig(dir+'/'+tag+'_boxplot_stripplot_5000_genes.pdf')
plt.close()


#------------------------------------------------------------------
# Generate barplot for tissue relationship
#------------------------------------------------------------------

decon2 = decon.merge(pheno[['tissue']], left_on = 'Mixture', right_on = 'Subject')
decon3 = decon2.merge(pheno[['Subtype']], left_on = 'Mixture', right_on = 'Subject')
decon2.set_index('tissue', inplace=True)
decon2.to_csv(dir+"/CIBERSORTx_ALS5000_with_tissue.csv")
decon2 = decon2.drop('Mixture', axis=1)
deconDF = decon2.T
deconDF.rename({'Lateral Motor Cortex':'Motor Cortex', 'Medial Motor Cortex':'Motor Cortex', 'Other Motor Cortex': 'Motor Cortex'}, axis=1, inplace= True)
celltypes1 = ['Glial_Progenitor','Excitatory','Neuronal_Progenitor','Inhibitory','Mural','Endothelial','Astrocyte', 'Microglia', 'Glyc', 'U1', 'Choroid']
deconDF = deconDF.loc[[True if i in celltypes1 else False for i in list(deconDF.index)]]
deconDF2 = deconDF.T
# Only specific subtypes
#subtypes1 = ['OX', 'TE', 'GLIA']
#deconDF2 = deconDF2.loc[[True if i in subtypes1 else False for i in list(deconDF2.index)]]
deconDF2 = deconDF2.reset_index()
deconDF2 = deconDF2[['tissue', 'Neuronal_Progenitor', 'Excitatory', 'Inhibitory', 'Glial_Progenitor', 'Astrocyte', 'Microglia', 'Endothelial', 'Mural', 'Choroid', 'U1']]
df_long = pd.melt(deconDF2, id_vars=['tissue'])
s1 = df_long.tissue.unique().tolist()
# Only specific subtypes
#s1_colors=["#b22222","#000080","#daa520"]
s1_colors=["#8571b2", "#F28C28"]
color_dict = dict(zip(s1, s1_colors))

# Plot box plot and stripplot
plt.figure(figsize=(15,10))
sns.set_theme(style="whitegrid")
sns.boxplot(x='variable', y='value', hue = 'tissue', linewidth=0.7, fliersize=2, palette = color_dict, data=df_long)
sns.stripplot(x='variable', y='value', hue = 'tissue', dodge=True, size=2, palette = color_dict, data=df_long)
plt.xticks(rotation = 30)
plt.ylim([0.00, 0.50])
plt.savefig(dir+'/'+tag+'_tissue_boxplot_stripplot_5000_genes.pdf')
plt.close()
