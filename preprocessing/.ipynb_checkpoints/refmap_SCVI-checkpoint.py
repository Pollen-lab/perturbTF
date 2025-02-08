#!/usr/bin/env python
# coding: utf-8

#Script to integrate data from all screens and reference map to Wang et al 2025 atlas

# In[3]:


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


# In[2]:


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


'''
# ## Li/Kreigstein multiome/RNA as reference

# ### Preprocessing/Subsetting original data

# In[4]:


directory='/wynton/group/pollen/jding/brainchromatin/Li/'
adata_ref = sc.read(os.path.join(directory,'adata.h5ad'))
adata_ref = adata_ref.raw.to_adata()
adata_ref.layers["counts"] = adata_ref.X.copy()
adata_ref.raw = adata_ref


# In[7]:


sc.pp.normalize_total(adata_ref,exclude_highly_expressed=True)
sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=2500, batch_key="Ident", subset= True)
#adata_ref.var['highly_variable']=[x for x in adata_ref.var['highly_variable'] if x in adata_query.var_names]
adata_ref.write(os.path.join(directory,'adata230k.h5ad'), compression='gzip')
'''

# ### Import query data

# In[61]:


directory='/wynton/group/pollen/jding/brainchromatin/Li/'
adata_ref = sc.read(os.path.join(directory,'ref.h5ad'), compression='gzip')


# In[9]:
adata1 = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/Cellbender/HD0.h5ad')
adata1.obs['timepoint'] = 'D0'
adata1.obs['species']='human'
adata2 = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/Cellbender/2D.h5ad')
adata2.obs['timepoint'] = 'D7'
adata2.obs['species']='human'
adata3 = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/Cellbender/HM2D.h5ad')
adata3.obs['timepoint'] = 'D7'
adata4 = sc.read('/wynton/group/pollen/jding/brainchromatin/HM2D/Cellbender/HM2D_2nd.h5ad')
adata4.obs['timepoint'] = 'D7'
#adata_query=anndata.AnnData.concatenate(adata1,adata2,adata3,adata4, join='outer', batch_categories=['HD0','HD7','HMD7','HM2ndD7']) 
adata_query=anndata.AnnData.concatenate(adata1,adata2,adata3,adata4, join='outer') 
#adata_query.obs['stage'] = [x.split('#')[0] for x in adata_query.obs['individual'].to_list()]
#adata_query.obs['sample_id'] = [x.split('-')[0] for x in adata_query.obs['individual'].to_list()]
#adata_query.obs['sex'] = [x.split('-')[1] for x in adata_query.obs['individual'].to_list()]


# In[11]:


adata_query.raw = adata_query
adata_query.layers['counts'] = adata_query.X
adata_ref.var_names_make_unique()
adata_query.var_names_make_unique()


# In[13]:


sc.pp.normalize_total(adata_query, target_sum=1e4,exclude_highly_expressed=True)
sc.pp.log1p(adata_query)
sc.pp.highly_variable_genes(adata_query, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata_query, max_value=10)
sc.pp.highly_variable_genes(adata_query, n_top_genes=2500, batch_key="batch_name", subset= False)


for x in adata_query.obs.columns.tolist():
    adata_query.obs[x]=adata_query.obs[x].astype(str)
for x in adata_query.var.columns.tolist():
    adata_query.var[x]=adata_query.var[x].astype(str)
adata_query.write(os.path.join(directory,'query_HM2Dall.h5ad'), compression='gzip')


# In[14]:


var = [x for x in adata_ref.var_names if x in adata_query.var_names]
adata_ref = adata_ref[:, var].copy()
adata_query = adata_query[:, var].copy()

'''
# In[15]:


adata_ref


# In[ ]:


scvi.model.SCVI.setup_anndata(adata_ref, batch_key="Ident", layer="counts")

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


# In[ ]:


adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)


# save the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
dir_path = os.path.join(directory,'scVI_vivo_model_HM2D/')
vae_ref.save(dir_path, overwrite=True)
'''

# In[ ]:


# load the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
dir_path = os.path.join(directory,'scVI_vivo_model_HM2D/')
vae_ref = scvi.model.SCVI.load(dir_path,adata_ref)


# In[ ]:


# both are valid
scvi.model.SCVI.prepare_query_anndata(adata_query, dir_path)
#scvi.model.SCVI.prepare_query_anndata(adata_query, vae_ref)


# In[ ]:

adata_query.obs['timepoint_batch_name'] = adata_query.obs['timepoint'].astype(str) + '-' + adata_query.obs['batch_name'].astype(str)
#adata_query.obs['Ident'] = adata_query.obs['batch_name']
adata_query.obs['Ident'] = adata_query.obs['timepoint_batch_name']


# In[ ]:


# both are valid
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    dir_path,
)
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    vae_ref,
)


# In[ ]:


vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()


# In[ ]:


sc.pp.neighbors(adata_query, use_rep="X_scVI")
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)


# In[ ]:


sc.pl.umap(adata_query, color=['timepoint','batch_name','individual','Ident'],frameon=False,ncols=1,save='scvi-query-all')


# ### Reference mapping from D40Human with SCANVI

# In[30]:s


adata_ref.obs["labels_scanvi"] = adata_ref.obs["type"].values


# In[31]:


directory='/wynton/group/pollen/jding/brainchromatin/Li/'
dir_path = os.path.join(directory,'scVI_vivo_model_HM2D/')
vae_ref = scvi.model.SCVI.load(dir_path, adata_ref)

'''
# In[32]:


# unlabeled category does not exist in adata.obs[labels_key]
# so all cells are treated as labeled
vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category="Unknown",
    labels_key="labels_scanvi",
)


# In[33]:


vae_ref_scan.train(max_epochs=20, n_samples_per_label=100)


# In[34]:


adata_ref.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)



# In[49]:


# save the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
dir_path_scan = os.path.join(directory,'scanvi_vivo_model_HM2D/')
vae_ref_scan.save(dir_path_scan, overwrite=True)
'''

# In[50]:


# load the reference model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
dir_path_scan = os.path.join(directory,'scanvi_vivo_model_HM2D/')
vae_ref_scan = scvi.model.SCANVI.load(dir_path_scan,adata_ref)


# In[51]:


# again a no-op in this tutorial, but good practice to use
scvi.model.SCANVI.prepare_query_anndata(adata_query, dir_path_scan)


# In[53]:


vae_q = scvi.model.SCANVI.load_query_data(adata_query,vae_ref_scan)


# In[54]:


vae_q.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)


# In[28]:


# save the query model
directory='/wynton/group/pollen/jding/brainchromatin/Li/'
dir_path_scan = os.path.join(directory,'scanvi_HM2Dall_model/')
vae_q.save(dir_path_scan, overwrite=True)


# In[55]:


adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
adata_query.obs["predictions"] = vae_q.predict()
#adata_query.write(os.path.join(directory,'query_HM2Dall.h5ad'), compression='gzip')


# In[56]:


sc.tl.leiden(adata_query,resolution=.8)
sc.pl.umap(adata_query,color=["predictions",'leiden'],frameon=False, legend_loc='on data',legend_fontsize ='xx-small')
adata_query.obs["leiden-.8"] = adata_query.obs["leiden"]


# In[57]:


sc.tl.leiden(adata_query,resolution=.3)
sc.pl.umap(adata_query,color=["predictions",'leiden'],frameon=False, legend_loc='on data',legend_fontsize ='xx-small')
adata_query.obs["leiden-.3"] = adata_query.obs["leiden"]


# In[30]:


sc.tl.leiden(adata_query,resolution=.2)
sc.pl.umap(adata_query,color=["predictions",'leiden'],frameon=False, legend_loc='on data',legend_fontsize ='xx-small')
adata_query.obs["leiden-.2"] = adata_query.obs["leiden"]


# In[58]:

sc.pl.umap(adata_query,color=["predictions",'leiden-.3','individual','Ident'],frameon=False,ncols=1,save='scanvi-query-all')


# In[38]:

viewgenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','BCAS1','HOPX','ERBB4','FOXG1','CALB1','CALB2','GAD1','GAD2',
            'GADD45G','LHX1','LHX5','LHX6','LHX8','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5',
            'NR2F1','NR2F2','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TH','DDC','SLC18A2','PAX6','MEIS2','ISL1','PENK','CRABP1',
            'ZIC1','ZIC4','EBF1','DLX5','TFAP2A','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','NEUORD2','NEUROG2','NEUROD1',
            'BCL11B','TLE4','FEZF2','SATB2','TBR1','SLC17A6','SLC17A7','RORB','AIF1','RELN','PECAM1','HBZ','TSHZ1','FOXP2','FOXP1']
viewgenes=[x for x in viewgenes if x in adata_query.var.index]
sc.pl.umap(adata_query,color=viewgenes,use_raw=False,save='genes-all')

import pandas as pd
adata_query.obs.to_csv(os.path.join(directory,'query_HM2Dall.obs.csv'))
df = pd.DataFrame(adata_query.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
df['cell_id'] = adata_query.obs.index
df.to_csv(os.path.join(directory,'query_HM2Dall.X_umap.csv'))
df = pd.DataFrame(adata_query.obsm['X_scANVI'])
df['cell_id'] = adata_query.obs.index
df.to_csv(os.path.join(directory,'query_HM2Dall.X_scANVI.csv'))

# Specify the filename to save the Scanpy object
filename = 'query_HM2Dall.pkl'

# Save the Scanpy object using pickle
with open(filename, 'wb') as f:
    pickle.dump(adata, f)


for x in adata_query.obs.columns.tolist():
    adata_query.obs[x]=adata_query.obs[x].astype(str)
for x in adata_query.var.columns.tolist():
    adata_query.var[x]=adata_query.var[x].astype(str)
adata_query.write(os.path.join(directory,'query_HM2Dall.h5ad'), compression='gzip')





