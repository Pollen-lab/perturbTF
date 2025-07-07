#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
#from scar import setup_anndata
import warnings
import anndata
import numpy as np
import os
#import bbknn
warnings.simplefilter("ignore")


# In[ ]:


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import anndata
import scvi
import scanpy as sc

sc.set_figure_params(figsize=(4, 4))
scvi.settings.seed = 94705


# ## 2D data as reference

# ### Preprocessing/Subsetting original data

# In[ ]:


adata_ref = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/D7-filtered_guides.IN.h5ad')


# In[ ]:


adata_ref


# In[ ]:


adata_ref.layers["counts"] = adata_ref.X.copy()
adata_ref.raw = adata_ref


# In[ ]:


sc.pp.normalize_total(adata_ref,exclude_highly_expressed=True)
sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=2500, batch_key="Ident", subset= True)
#adata_ref.var['highly_variable']=[x for x in adata_ref.var['highly_variable'] if x in adata_query.var_names]


# In[8]:


sc.pl.umap(
    adata_ref,
    color=[ "individual",'subtype'],
    frameon=False,
)


# In[61]:


adata_ref.obs['subtype'].value_counts()



# ### Import query data

# In[3]:


#subset data for 2D IN lineage
import scanpy as sc 
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
adata_query = sc.read(os.path.join(directory,'query_LMO1.h5ad'), compression='gzip')
adata_query.obs['perturbation'] = ['NT' if x in ['WT','non-targeting'] else 'Perturbed' for x in adata_query.obs['gene_NKS']]
adata_query.obs['supervised_name'] = adata_query.obs['predictions']
adata_query = adata_query[adata_query.obs['predictions'].isin(['IN_dLGE-CGE'])]
adata_query = adata_query.raw.to_adata()

adata2 = sc.read(os.path.join(directory,'HM2D_IN_query.h5ad'), compression='gzip')
#adata2 = adata2[adata2.obs['species'] == 'human']
adata_query = anndata.AnnData.concatenate(adata_query, adata2, join='outer', batch_categories=['LMO1','HM2D']) 
adata_query.obs['batch'] = adata_query.obs['batch'].astype(str)

print(adata_query)
print(adata_query.obs['batch'])

# In[ ]:


#adata_query = adata_query[adata_query.obs['species'] =='macaque']


# In[9]:


adata_query.raw = adata_query.copy()
adata_query.layers['counts'] = adata_query.X


# In[10]:


adata_ref.var_names_make_unique()
adata_query.var_names_make_unique()

print(adata_query)

# In[11]:


sc.pp.normalize_total(adata_query, target_sum=1e4,exclude_highly_expressed=True)
sc.pp.log1p(adata_query)
#sc.pp.highly_variable_genes(adata_query, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.pp.scale(adata_query, max_value=10)
#sc.pp.highly_variable_genes(adata_query, n_top_genes=2500, batch_key="individual", subset= False)


# In[12]:


var = [x for x in adata_ref.var_names if x in adata_query.var_names]
adata_ref = adata_ref[:, var].copy()
adata_query = adata_query[:, var].copy()


# In[13]:


scvi.model.SCVI.setup_anndata(adata_ref, batch_key="individual", layer="counts")

arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)

vae_ref = scvi.model.SCVI(
    adata_ref,
    **arches_params
)
vae_ref.train()


# In[14]:


adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)


# In[15]:


sc.pl.umap(
    adata_ref,
    color=[ "individual",'subtype'],
    frameon=False,
)


# In[16]:


# save the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
directory='/wynton/scratch/jding/brainchromatin/Li/'
dir_path = os.path.join(directory,'scVI_model_IN/')
vae_ref.save(dir_path, overwrite=True)


# In[17]:


# load the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
directory='/wynton/scratch/jding/brainchromatin/Li/'
dir_path = os.path.join(directory,'scVI_model_IN/')
vae_ref = scvi.model.SCVI.load(dir_path,adata_ref)


# In[18]:


# both are valid
scvi.model.SCVI.prepare_query_anndata(adata_query, dir_path)
#scvi.model.SCVI.prepare_query_anndata(adata_query, vae_ref)


# In[19]:


#adata_query.obs['individual'] = adata_query.obs['batch']


# In[20]:


# both are valid
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    dir_path,
)
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    vae_ref,
)


# In[21]:


vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()


# In[22]:


sc.pp.neighbors(adata_query, use_rep="X_scVI")
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)


# In[23]:


sc.pl.umap(
    adata_query,
    color=[ "individual"],
    frameon=False,
)



# ### Reference mapping from 2D with SCANVI

# In[24]:


adata_ref.obs["labels_scanvi"] = adata_ref.obs["subtype"].values


# In[25]:


directory='/wynton/group/pollen/jding/brainchromatin/Li/'
directory='/wynton/scratch/jding/brainchromatin/Li/'
dir_path = os.path.join(directory,'scVI_model_IN/')
vae_ref = scvi.model.SCVI.load(dir_path, adata_ref)


# In[26]:


# unlabeled category does not exist in adata.obs[labels_key]
# so all cells are treated as labeled
vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",
    labels_key="labels_scanvi",
)


# In[27]:


vae_ref_scan.train(max_epochs=20, n_samples_per_label=100)


# In[28]:


adata_ref.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)


# In[30]:


# save the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
directory='/wynton/scratch/jding/brainchromatin/Li/'
dir_path_scan = os.path.join(directory,'scanvi_model_IN/')
vae_ref_scan.save(dir_path_scan, overwrite=True)


# In[31]:


# load the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
directory='/wynton/scratch/jding/brainchromatin/Li/'
dir_path_scan = os.path.join(directory,'scanvi_model_IN/')
vae_ref_scan = scvi.model.SCANVI.load(dir_path_scan,adata_ref)


# In[32]:


# again a no-op in this tutorial, but good practice to use
scvi.model.SCANVI.prepare_query_anndata(adata_query, dir_path_scan)


# In[33]:


vae_q = scvi.model.SCANVI.load_query_data(adata_query,vae_ref_scan)


# In[34]:


vae_q.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)


# In[35]:


# save the query model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
directory='/wynton/scratch/jding/brainchromatin/Li/'
dir_path_scan = os.path.join(directory,'scanvi_model_IN_query/')
vae_q.save(dir_path_scan, overwrite=True)


# In[36]:


adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
adata_query.obs["predictions"] = vae_q.predict()


# In[37]:


sc.pl.umap(adata_query,color=["predictions",'individual','gene_NKS'],frameon=False,ncols=2, save = 'IN')
sc.pl.umap(adata_query,color=["predictions",'individual','gene_NKS'],frameon=False,ncols=2, save = 'IN')


clonal_cluster_count = adata_query.obs.groupby(["gene_NKS","predictions"]).size().reset_index().pivot(
    index="gene_NKS", columns="predictions", values=0
)
clonal_cluster_count = clonal_cluster_count[clonal_cluster_count.columns[clonal_cluster_count.columns != "NA"]]
clonal_cluster_count = clonal_cluster_count[clonal_cluster_count.sum(axis=1) > 5]
clonal_cluster_count = (clonal_cluster_count.T / clonal_cluster_count.sum(axis=1)).T

fig, ax = plt.subplots(figsize=(12, 6))
ax = clonal_cluster_count.plot(kind="bar", stacked=True, width=0.9, edgecolor="black",
                               #color=dict(zip(clones.obs.gene_NKS.cat.categories, clones.uns["perturbation_colors"])),
                               ax=ax, linewidth=0.5)
ax.legend(loc=(1.1, 0))
ax.set_xlabel("")
fig.savefig("INsubtype.pdf")


import pandas as pd
adata_query.obs.to_csv(os.path.join(directory,'LMO1_IN_query.obs.csv'))
df = pd.DataFrame(adata_query.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
df['cell_id'] = adata_query.obs.index
df.to_csv(os.path.join(directory,'LMO1_IN_query.X_umap.csv'))
df = pd.DataFrame(adata_query.obsm['X_scANVI'])
df['cell_id'] = adata_query.obs.index
df.to_csv(os.path.join(directory,'LMO1_IN_query.X_scANVI.csv'))


# In[39]:
for x in adata_query.obs.columns.tolist():
    adata_query.obs[x]=adata_query.obs[x].astype(str)
for x in adata_query.var.columns.tolist():
    adata_query.var[x]=adata_query.var[x].astype(str)

for x in adata_query.obs.columns.tolist():
    adata_query.obs[x]=adata_query.obs[x].astype(str)
for x in adata_query.var.columns.tolist():
    adata_query.var[x]=adata_query.var[x].astype(str)
    
adata_query.obs['mt-LMO1'] = adata_query.obs['mt-LMO1'].astype(str)

adata_query.write('/wynton/group/pollen/jding/brainchromatin/Li/LMO1_IN_query.h5ad')


