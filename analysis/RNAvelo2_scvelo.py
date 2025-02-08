#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import anndata
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import re
import sklearn
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from collections import Counter
import random
import seaborn
import sys
import shutil
import scvelo as scv
import tqdm
from ipywidgets import IntProgress

# In[3]:


sc.settings.figdir='figures/scVelo/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True
scv.settings.figdir='figures/scVelo/'
scv.set_figure_params(style='scvelo',format='pdf')


import os 
directory='/wynton/group/pollen/jding/brainchromatin/perturb/velocyto/'
namesOfInterest=['2D_L2', '2D_L3', '2D_L4']


# In[17]:

import anndata
import scipy
import pandas as pd

#Load annotated cellbended output
adata = sc.read('/wynton/group/pollen/jding/brainchromatin/perturb/Cellbender/2D.h5ad')

#load VeloCyto output
adatas=[]
for f in namesOfInterest:
    print(f)
    spliced = scipy.io.mmread(os.path.join(directory,f,'Velocyto/filtered/spliced.mtx')).tocsr()
    unspliced = scipy.io.mmread(os.path.join(directory,f,'Velocyto/filtered/unspliced.mtx')).tocsr()
    ambiguous = scipy.io.mmread(os.path.join(directory,f,'Velocyto/filtered/ambiguous.mtx')).tocsr()
    obs=pd.read_csv(os.path.join(directory,f,'Velocyto/filtered/barcodes.tsv'),index_col=0,header=None)
    var=pd.read_csv(os.path.join(directory,f,'Velocyto/filtered/features.tsv'),index_col=1,sep='\t',header=None)
    ad = anndata.AnnData(obs= obs, var= var, layers={'spliced': spliced.T,'unspliced': unspliced.T, 'ambiguous': ambiguous.T})
    ad.obs['batch_name']=f
    ad.obs['barcode']=ad.obs_names
    ad.var.index.name = "GeneID"
    ad.var.columns =['Symbol', 'Type']
    ad.obs.index.name = "Barcode"
    ad.obs_names_make_unique()
    ad.var_names_make_unique()
    adatas.append(ad)

ad=anndata.AnnData.concatenate(*adatas, join='outer', batch_categories=namesOfInterest) 
ad.obs_names = ad.obs['barcode']+'-1-'+ad.obs['batch_name']
ad

#subset cells/genes existing in both dataset
ad.obs_names_make_unique()
ad.var_names_make_unique()    
ad = ad[ad.obs_names.isin(adata.obs_names)]
adata = adata[adata.obs_names.isin(ad.obs_names)]
ad = ad[:,ad.var_names.isin(adata.var_names)]
adata = adata[:,adata.var_names.isin(ad.var_names)]
ad = ad[adata.obs_names,:]
ad = ad[:,adata.var_names]
print(ad)

#transfer splic/unsplice matrix to cellbended object
adata.layers['spliced'] = ad.layers['spliced']
adata.layers['unspliced'] = ad.layers['unspliced']
adata.layers['ambiguous'] = ad.layers['ambiguous']

print(adata)
#subset to non-targeting ctrls
adata = adata[adata.obs['num_gene_IDs'].isin([1,2])]
adata.obs['Gene_target'] = adata.obs['gene_IDs']
adata.obs['Gene_target'] =[x.replace('non-targeting,','') for x in adata.obs['Gene_target']] 
adata.obs['Gene_target'] =[x.replace(',non-targeting','') for x in adata.obs['Gene_target']] 
adata = adata[adata.obs['Gene_target'] =='non-targeting']
print(adata)


#adata=sc.read(os.path.join(directory, '2D_scVelo.h5ad'))
scv.pp.filter_genes(adata, min_shared_counts=20)

#calculate variable genes
sc.pp.highly_variable_genes(adata,flavor='seurat_v3',layer='unspliced',n_top_genes=5000)
adata.var['highly_variable_rank_u']=adata.var['highly_variable_rank']
sc.pp.highly_variable_genes(adata,flavor='seurat_v3',layer='spliced',n_top_genes=5000)
adata.var['highly_variable_rank_s']=adata.var['highly_variable_rank']
adata.var['velocity_genes']=adata.var.loc[:,['highly_variable_rank_s','highly_variable_rank_u']].mean(1,skipna=False).rank()<3000


# In[10]:


scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)
scv.pp.remove_duplicate_cells(adata)
scv.pp.moments(adata, n_neighbors=20)
scv.tl.recover_dynamics(adata,var_names='velocity_genes',use_raw=False, n_jobs = 16)
scv.tl.velocity(adata,mode='dynamical',filter_genes=False)
scv.tl.velocity_graph(adata,mode_neighbors='connectivities',vkey='velocity',approx=False, n_jobs = 16)
scv.tl.latent_time(adata)

adata.write('/wynton/scratch/jding/brainchromatin/HM2D/scVelo-D7NT.h5ad')



# In[ ]:
scv.pl.velocity_embedding_stream(adata, basis='umap', save="leiden-ondata"+gene, color = "supervised_name")
scv.pl.velocity_embedding_stream(adata, basis='umap', save="leiden-margin"+gene, color = "leiden", legend_loc='right margin')
scv.pl.velocity_embedding_stream(adata, basis='umap', save="individual"+gene, color = "individual")
scv.pl.velocity_embedding_stream(adata, basis='umap', save="class"+gene, color = "class")

