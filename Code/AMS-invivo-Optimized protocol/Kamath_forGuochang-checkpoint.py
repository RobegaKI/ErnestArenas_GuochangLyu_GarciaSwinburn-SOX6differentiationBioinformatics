#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import scanpy as sc
import anndata as ad
import loompy


# In[2]:


kamath = '/datb/ea/Carmen/atlas_proj/adult/loom_files/kamath/kamath_ctrl.loom'


# In[3]:


with loompy.connect(kamath, 'r') as ds:
    print(ds.ca.keys())
    print(ds.ra.keys())


# In[4]:


adata_kamath = sc.read_loom(kamath,
    sparse=True,
    X_name="counts",
    obs_names="Cell_ID",
    var_names="Gene")


# In[6]:


sc.pl.embedding(adata_kamath, basis = 'X_umap', color = 'Cell_Type', cmap = 'inferno_r', frameon = False)


# In[7]:


kamath_CT = '/datb/ea/Carmen/atlas_proj/adult/20250128_version/celltypist_labeltransfer/20250115_Kamath-nonintegrated-CT.loom'


# In[8]:


adata_kamath_CT = sc.read_loom(kamath_CT,
    sparse=True,
    X_name="counts",
    obs_names="Cell_ID",
    var_names="Gene")


# In[9]:


adata_kamath_CT


# In[12]:


sc.pl.embedding(adata_kamath_CT, basis = 'UMAP_paper', color = 'majority_voting', cmap = 'inferno_r', frameon = False)


# In[14]:


print(adata_kamath.shape[0])
print(adata_kamath_CT.shape[0])


# In[21]:


adata_kamath_CT.obs_names = adata_kamath_CT.obs['CellID']
adata_kamath.obs_names = adata_kamath.obs['CellID']


# In[26]:


majority_ = dict(zip(adata_kamath_CT.obs['CellID'].copy(), adata_kamath_CT.obs['majority_voting'].copy()))
conf_score = dict(zip(adata_kamath_CT.obs['CellID'].copy(), adata_kamath_CT.obs['conf_score'].copy()))
pred_labels = dict(zip(adata_kamath_CT.obs['CellID'].copy(), adata_kamath_CT.obs['predicted_labels'].copy()))


# In[33]:


majority_full = dict(zip(adata_kamath.obs['CellID'].copy(), np.repeat('NotTrained', adata_kamath.shape[0])))
conf_score_full = dict(zip(adata_kamath.obs['CellID'].copy(), np.repeat(np.nan, adata_kamath.shape[0])))
pred_labels_full = dict(zip(adata_kamath.obs['CellID'].copy(), np.repeat('NotTrained', adata_kamath.shape[0])))


# In[34]:


for cell in majority_.keys():
    majority_full[cell] = majority_[cell]
    conf_score_full[cell] = conf_score[cell]
    pred_labels_full[cell] = pred_labels[cell]


# In[35]:


adata_kamath.obs['majority_voting'] = majority_full
adata_kamath.obs['conf_score'] = conf_score_full
adata_kamath.obs['predicted_labels'] = pred_labels_full


# In[37]:


sc.pl.embedding(adata_kamath, basis = 'X_umap', color = ['majority_voting', 'predicted_labels', 'conf_score'],
                cmap = 'inferno_r', wspace = 0.3, frameon = False)


# In[38]:


get_ipython().system("mkdir '/proj/user/carmen/forGuochang/'")
adata_kamath.write('/proj/user/carmen/forGuochang/Kamath-CellTypist.h5ad')


# In[ ]:




