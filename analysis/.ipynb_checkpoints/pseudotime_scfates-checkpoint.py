#!/usr/bin/env python
# coding: utf-8

#script to compute pseudotime and branches in non-targeting cells on differentiation D7 

# ### Loading modules and settings

# In[1]:


import os, sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
import scFates as scf
import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import numpy as np
import random
import pandas as pd
import scanpy as sc
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
#sc.settings.logfile = sys.stdout
## fix palantir breaking down some plots
import seaborn 
seaborn.reset_orig()

sc.set_figure_params()
scf.set_figure_pubready()


# ## Preprocessing pipeline from Palantir

# ### Load, normalize and log-transform count data


# In[42]:

adata = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/D7-filtered_guides.h5ad', compression='gzip')
adata.obs['class'] = adata.obs['class'].str.replace(' ', '')
adata.obs['Gene_target'] = adata.obs['gene_NKS'] 
adata = adata[adata.obs['Gene_target'] == 'non-targeting']
adata.var_names_make_unique()
adata.obs_names_make_unique()

sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500)
sc.pp.pca(adata)


# ### Perform PCA on highly variable genes

# In[45]:


pca_projections = pd.DataFrame(adata.obsm["X_pca"],index=adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# ### Generate embedding from the multiscale diffusion space


# generate neighbor draph in multiscale diffusion space
adata.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adata,n_neighbors=30,use_rep="X_palantir")


# In[65]:


# draw ForceAtlas2 embedding using 2 first PCs as initial positions
adata.obsm["X_pca2d"]=adata.obsm["X_pca"][:,:2]
sc.tl.draw_graph(adata,init_pos='X_pca2d')



# ## Tree learning with SimplePPT

# In[68]:

#increase ppt_lambda helps make the trajectory smooth
scf.tl.tree(adata,method="ppt",Nodes=50,use_rep="palantir",
            device="cpu",seed=1,ppt_lambda=300,ppt_sigma=0.025,ppt_nsteps=600)


# ### projecting results onto ForceAtlas2 embedding

# In[69]:


scf.pl.graph(adata)


# ### Selecting a root and computing pseudotime

# In[70]:


scf.tl.root(adata,18)
scf.tl.pseudotime(adata,n_jobs=20,n_map=100,seed=42)


# ## Representing the trajectory and tree
# 
# ### on top of existing embedding

# In[13]:


scf.pl.trajectory(adata)



# In[14]:


sc.pl.draw_graph(adata,color=["seg","milestones"])
scf.tl.rename_milestones(adata,["RG_DIV","RG/Astro", "IN", 'EN'])
# we change the color of the root milestone for better visualisations
adata.uns["milestones_colors"][3]="#17bece" 




# ### as a dendrogram representation

# In[17]:


scf.tl.dendrogram(adata)


# In[18]:


scf.pl.dendrogram(adata,color="seg")


# In[19]:


sc.set_figure_params(figsize=(1.5,4),frameon=False,dpi_save=300)
scf.pl.dendrogram(adata,color="t",show_info=False,save="tree1.pdf",cmap="viridis")
scf.pl.dendrogram(adata,color="milestones",legend_loc="on data",color_milestones=True,legend_fontoutline=True,save="tree2.pdf")


# ## Test and fit features associated with the tree
# 
# Let's find out which genes are significantly changing along the tree.

# In[20]:


scf.tl.test_association(adata,n_jobs=20)
scf.tl.fit(adata,n_jobs=20)

# In[21]:


sc.set_figure_params()
scf.pl.test_association(adata)
plt.savefig("figures/C.pdf",dpi=300)

adata.write('/wynton/group/pollen/jding/brainchromatin/HM2D/D7-filtered_guides.h5ad')


#plot TF expression along pseudotime
genes= ["NR4A2", "POU2F1", "SOX5", "VEZF1", "EMX1", "BHLHE22", "ZBTB20", "SOX2" , "TCF7L1", "TFAP2C" ,"MEF2C", "POU3F2" , "ARX", "TCF3", "ZBTB18" , "SOX6" ,
  "SOX9" , "NR2E1", "CTCF", "NEUROD6", "JUN", "ZNF148", "NFIB" , "ZNF219", "NEUROD2" , "ZNF281", "NFIA", "E2F1", "ZNF441", "FOS", "TCF12", "POU3F1", "TBR1", "KLF11", "KLF10",
  "SATB2", "KAT7", "ATF7", "RFX2", "KLF3", "MEIS2", "PHF21A", "DEAF1",  "ASCL1"]
genes = [x for x in genes if x in adata.var_names]
len(genes)

scf.tl.fit(adata,features=genes, n_jobs=4)

g=scf.pl.trends(adata, 
                 features=genes,
                 highlight_features=[], n_features=0,
                 plot_emb=False,ordering="max",return_genes=True)

sc.set_figure_params()
scf.pl.matrix(adata,g,norm="minmax",cmap="RdBu_r",#annot_var=True,
              colorbar=True,save="_perturb.pdf")

