#!/usr/bin/env python
# coding: utf-8

# This notebook is a tutorial for computing SEACell metacells, visualizing results and computing evaluation metrics

# # Imports

# In[8]:


import numpy as np
import pandas as pd
import scanpy as sc
import os


# In[2]:


import SEACells


# In[3]:


import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# In[4]:


# Some plotting aesthetics
#get_ipython().run_line_magic('matplotlib', 'inline')

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100


# # Load Data
# 
# We recommend the use of scanpy Anndata objects as the preferred mode of loading and filtering data.
# 
# A sample datset is available for download with the instructions listed below. This is a filtered, unnormalized counts of single-nuclear RNA-seq dataset of CD34+ sorted bone marrow cells to profile human hematopoiesis [Dataset ref TBD].
# 
# Uncomment the following lines to download the sample dataset in a Unix-based system. For non-UNIX systems, download the files using the URL

# In[5]:


#load data 
import os
work_dir = '/wynton/group/pollen/jding/brainchromatin/Li/scenicplus'

if not os.path.exists(os.path.join(work_dir, 'seacells')):
    os.makedirs(os.path.join(work_dir, 'seacells'))
data_dir = os.path.join(work_dir, 'seacells')

ad = sc.read(os.path.join(work_dir, 'scATAC/cistopic_obj.h5ad'))

# In[6]:


print(ad)


# In[81]:


directory = os.path.expanduser('/wynton/group/pollen/jding/brainchromatin/Li/')
adata = sc.read(os.path.join(directory,'adata.h5ad'))
adata.obs.index =  adata.obs['gex_barcode'].astype(str) +'___'+ adata.obs['batch_name'].astype(str) 



# In[ ]:


adata = adata[adata.obs_names.isin(ad.obs_names)]
ad = ad[adata.obs_names,:]
print(ad)
print(adata)


# In[85]:


# Plot cell-types for reference
for x in adata.obs.columns:
    ad.obs[x] = adata.obs[x]
ad.obsm['X_umap'] = adata.obsm['X_umap'] 
sc.pl.umap(ad,  color='subclass', frameon=False,legend_loc='on data',legend_fontsize=8,save='subclass')
sc.pl.umap(ad,  color='type', frameon=False,legend_loc='on data',legend_fontsize=8,save='type')

# # Running SEACells

# As a rule of thumb, we recommended choosing one metacell for every 75 single-cells. Since this dataset contains ~7k cells, we choose 90 metacells.
# 
# <b>Note 1: </b> Running SEACells modifies the input Anndata object and adds the SEACell metacell assignments to the `obs` dataframe in the anndata object.
# <b>Note 2: </b> This analysis takes approxmiately 5 minutes

# In[11]:


## User defined parameters

## Core parameters 
n_SEACells = int(np.floor(ad.obs.shape[0] / 75))
build_kernel_on = 'X_topic' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells
use_gpu = True

# In[12]:


model = SEACells.core.SEACells(ad, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5,
                  use_gpu = use_gpu)


# In[13]:


model.construct_kernel_matrix()
#sp.sparse.save_npz(os.path.join(directory,'kernel_matrix.npz'), model.kernel_matrix)
#import scipy as sp
#model.kernel_matrix = sp.sparse.load_npz(os.path.join(directory,'kernel_matrix.npz'))

#import pickle
#model=pickle.load(open(data_dir + "/model.pkl",'rb'))

# In[15]:


# Initialize archetypes
model.initialize_archetypes()


# In[16]:


# Plot the initilization to ensure they are spread across phenotypic space
SEACells.plot.plot_initialization(ad, model,save_as='initilization.png')


# In[17]:


model.fit(min_iter=10, max_iter=50)

# In[18]:


# You can force the model to run additional iterations step-wise using the .step() function
#print(f'Ran for {len(model.RSS_iters)} iterations')
#for _ in range(5):
#    model.step()
#print(f'Ran for {len(model.RSS_iters)} iterations')


# # Accessing results

# ## Model Convergence

# In[19]:


# Check for convergence 
model.plot_convergence()


# ## SEACell Hard Assignments
# 
# These can be accessed as folows:
# - in the modified anndata object in `.obs['SEAell']` 
# - from the model using `.get_hard_assignments()` 
# 

# In[20]:


ad.obs[['SEACell']].head()


#In[21]:


model.get_hard_assignments().head()
ad.write(data_dir + '/seacells.h5ad')


ad.layers['raw'] = ad.X
SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')
SEACell_ad.write(data_dir + '/atac_meta.h5ad')

import pickle
with open(data_dir + "/model.pkl", "wb") as f:
            pickle.dump(model, f)
    

# In[25]:
adata.obs['SEACell'] = ad.obs['SEACell']
rna_meta_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
rna_meta_ad.write(os.path.join(data_dir, 'rna_meta_ad.h5ad'))

rna_meta_ad.obs['celltype'] = adata.obs.groupby('SEACell').apply(lambda x: pd.Series(x['type']).mode())#%%
rna_meta_ad.obs['subclass'] = adata.obs.groupby('SEACell').apply(lambda x: pd.Series(x['subclass']).mode())#%%
SEACell_ad.obs['celltype'] = ad.obs.groupby('SEACell').apply(lambda x: pd.Series(x['type']).mode())#%%
SEACell_ad.obs['subclass'] = ad.obs.groupby('SEACell').apply(lambda x: pd.Series(x['subclass']).mode())#%%
rna_meta_ad.write(data_dir + '/rna_meta.h5ad')
SEACell_ad.write(data_dir + '/atac_meta.h5ad')


# ## Normalization

# Normalization of metacell data can be performed using the `sc.pp.normalize_total` and `sc.pp.log1p` functions

# # Evaluating Results
# 
# We provide several methods for evaluating SEACell assignments:

# ## Visualizing Results
# 
# Metacells also implements methods for visualizing the results of the Metacells algorithm 
#     <ul> 
#         <li>```.plot_2D()``` provides an interface for viewing metacell assignments on any 2-dimensional embedding in ad.obsm. Plots can also be coloured by metacell assignment.
#         <li>```.plot_SEACell_sizes()``` can be used to view the distribution of number of cells assigned to each metacell
#     </ul>
#     
#             

# In[30]:


SEACells.plot.plot_2D(ad, key='X_umap', colour_metacells=False,save_as='umap_seacells.png')

# In[31]:


SEACells.plot.plot_2D(ad, key='X_umap', colour_metacells=True,save_as='umap_seacells_color.png')


# In[32]:


SEACells.plot.plot_SEACell_sizes(ad, bins=5,save_as='dis_metacells.png')


# ## Quantifying Results
# 
# SEACells also implements methods for visualizing the results of the SEACells algorithm 
#     <ul> 
#         <li>```.compute_celltype_purity(ad, col_name)``` computes the purity of different celltype labels within a SEACell metacell. Typically, col_name='celltype' or similar. Returns a pd.DataFrame of length n_SEACells.
#         <li>```.compactness(ad, low_dim_embedding)``` computes the per-SEAcell variance in diffusion components. ```low_dim_embedding``` is a string specifying the low dimensional embedding with which diffusion components are calculated, typically 'X_pca' for RNA or 'X_svd' for ATAC. Lower values of compactness suggest more compact/lower variance metacells.
#         <li>```separation(ad, low_dim_embedding,nth_nbr=1,cluster=None)``` computes the diffusion distance between a SEACell and its ```nth_nbr```. As before, ```low_dim_embedding``` is a string specifying the low dimensional embedding with which diffusion components are calculated, typically 'X_pca' for RNA or 'X_svd' for ATAC. If ```cluster``` is provided as a string, e.g. 'celltype', nearest neighbors are restricted to have the same celltype value.  Higher values of separation suggest better distinction between metacells.
#     </ul>
#     
# 

# In[33]:


SEACell_purity = SEACells.evaluate.compute_celltype_purity(ad, 'type')

plt.figure(figsize=(4,4))
sns.boxplot(data=SEACell_purity, y='type_purity')
plt.title('type Purity')
sns.despine()
plt.savefig("purity.png") 
plt.show()
plt.close()

SEACell_purity.head()


# In[34]:


compactness = SEACells.evaluate.compactness(ad, 'X_svd')

plt.figure(figsize=(4,4))
sns.boxplot(data=compactness, y='compactness')
plt.title('Compactness')
sns.despine()
plt.savefig("Compactness.png") 
plt.show()
plt.close()

compactness.head()


# In[35]:


separation = SEACells.evaluate.separation(ad, 'X_svd',nth_nbr=1)

plt.figure(figsize=(4,4))
sns.boxplot(data=separation, y='separation')
plt.title('Separation')
sns.despine()
plt.savefig("Separation.png") 
plt.show()
plt.close()

separation.head()


# In[36]:


ad.write(data_dir + '/seacells.h5ad')

