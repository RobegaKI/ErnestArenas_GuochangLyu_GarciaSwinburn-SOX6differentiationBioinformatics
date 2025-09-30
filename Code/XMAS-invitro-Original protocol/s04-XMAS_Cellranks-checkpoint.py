#!/usr/bin/env python
# coding: utf-8

# # 0. Preprocessing

# In[1]:


import cellrank as cr
import palantir


# In[2]:


import anndata
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt


# In[3]:


import os

# Set up global parameters
OS_path = "/Users/gclyu07/Desktop/XMAS_analysis/"
OS_path_datasets = os.path.join(OS_path, "dataset/")
OS_path_inputs = os.path.join(OS_path, "inputs/")
OS_path_outputs = os.path.join(OS_path, "outputs/")

seed = 1121


# In[4]:


os.getcwd()


# In[6]:


adata = sc.read_h5ad("./outputs/s04_XMAS_velocyto_spliced_unspliced_corrected.h5ad")


# In[7]:


adata


# In[8]:


adata.layers['spliced']


# In[9]:


adata.layers['unspliced']

cell_list = {
    0:'DA',
    1:'NbM',
    2:'Gaba',
    3:'GabaNb',
    4:'ProgFP',
    5:'Rgl1',
    6:'Rgl2',
    7:'Rgl3',}
adata.obs['bimod_oc'] = (
adata.obs['bimod_oc']
.astype('category')
.map(cell_list)
)cell_list = {
0:'Rgl1',
1:'Rgl2',
2:'Rgl3 ECM organization',
3:'Rgl3 Purine metabolism',
4:'Rgl3 Fate commitment',
5:'ProgFP S',
6:'ProgFP G2M',
7:'ProgFP Development',
8:'NbM',
9:'GabaNb',
10:'Gaba',
11:'DA Neuron projection',
12:'DA Synaptic modulation',
13:'DA Synaptic assembly',
14:'DA Neurotransmitter release'}
adata.obs['bimod_oc_des'] = (
adata.obs['bimod_oc_des']
.astype('category')
.map(cell_list)
)
# In[41]:


adata.write_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/s04_XMAS_velocyto_spliced_unspliced_corrected.h5ad")


# In[10]:


oi_colors={ "D11": "#F0AD4E", "D16": "#D9534F", "D28": "#428BCA", "D42": "#9933CC", "D56":"#66CCCC" }

adata.obs['bimod_oc_des'] = adata.obs['bimod_oc_des'].cat.rename_categories({'DA Postsynaptic transmission': 'DA Neurotransmitter transmission',})adata.obs['seurat_clusters_BiMod_OC'] = adata.obs['seurat_clusters_BiMod_OC'].astype(str)adata.obsm['umap']=adata.obsm['X_umap_pca_RNA_0.3']
# In[11]:


sc.set_figure_params(figsize=(4, 3))
sc.pl.umap(adata, color='orig.ident')
sc.pl.umap(adata, color='bimod_oc_des')
sc.pl.umap(adata, color='bimod_oc')


# In[45]:


sc.set_figure_params(figsize=(5, 4))
sc.pl.umap(adata, color=['orig.ident'], add_outline=True, palette=oi_colors,
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Orignal identities'])

adata.obsm['umap']=adata.obsm['X_umap_BiMod_0.3']
# In[41]:


sc.set_figure_params(figsize=(5, 4))
sc.pl.umap(adata, color=['orig.ident'], add_outline=True, palette=oi_colors,
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Orignal identities'])

adata.obsm['umap']=adata.obsm['X_umap_ATAC_0.3']
# In[43]:


sc.set_figure_params(figsize=(5, 4))
sc.pl.umap(adata, color=['orig.ident'], add_outline=True, palette=oi_colors,
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Orignal identities'])


# ___

# # 1. Palantir preprocessing

# In[12]:


adata = adata[(~adata.obs['bimod_oc_des'].isin(['Rgl1', 'Rgl2', 'Rgl3 ECM organization', 'Rgl3 Purine metabolism', 'Rgl3 Fate commitment',
                                               'ProgFP S', 'ProgFP G2M', 'Gaba', 'GabaNb']))]


# In[13]:


adata = adata[~(adata.obs['bimod_oc'].isin(['NbM','ProgFP']) & (adata.obs['orig.ident'].isin(['D28','D42','D56'])))]
adata = adata[~(adata.obs['bimod_oc'].isin(['ProgFP']) & (adata.obs['orig.ident'].isin(['D16'])))]


# In[14]:


sc.pp.neighbors(adata, n_neighbors=40, n_pcs=15, use_rep='X_pca', method='gauss')


# In[15]:


sc.tl.umap(adata)


# In[16]:


sc.set_figure_params(figsize=(5, 4))
sc.pl.umap(adata, color=['orig.ident','bimod_oc_des'], add_outline=True, legend_loc='on data',
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Orignal identities','Cell types'])


# In[17]:


dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)


# In[18]:


ms_data = palantir.utils.determine_multiscale_space(adata)


# In[19]:


imputed_X = palantir.utils.run_magic_imputation(adata)


# In[20]:


palantir.plot.plot_diffusion_components(adata)


# In[21]:


Tn_mask = np.isin(adata.obs['bimod_oc_des'], ['ProgFP Development'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask, 1])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata[[root_id]].obs_names


# In[22]:


palantir.plot.highlight_cells_on_umap(adata, pd.Series(["ProgFP Development"],
                index=["CCAGGATGTCAAAGGG-1"]))


# In[23]:


Tn_mask = np.isin(adata.obs['bimod_oc_des'], ['DA Synaptic modulation'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask,1])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata[[root_id]].obs_names


# In[24]:


start_cell = 'CCAGGATGTCAAAGGG-1'
terminal_states = pd.Series(["DA Neuron projection","DA Synaptic assembly",
                             "DA Neurotransmitter transmission","DA Synaptic modulation"],
                index=["CTTGTAAAGTAAAGGT-3", "TATGGATGTCACGAAC-5", "CCTTACTCATTATGAC-4", "ATTAGCGGTGCTGGTG-5"])


# In[25]:


palantir.plot.highlight_cells_on_umap(adata, terminal_states)


# In[26]:


pr_res = palantir.core.run_palantir(
    adata, start_cell, num_waypoints=500, terminal_states=terminal_states.index
)
pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]


# In[27]:


palantir.plot.plot_palantir_results(adata, s = 2.5)


# In[28]:


sc.pl.embedding(adata, basis="umap", color=["orig.ident","palantir_pseudotime","bimod_oc_des"], cmap="viridis_r")


# # 2. Kernels combination

# In[29]:


from cellrank.kernels import PseudotimeKernel
pk = PseudotimeKernel(adata, time_key="palantir_pseudotime")


# In[30]:


pk.compute_transition_matrix()


# In[31]:


scv.set_figure_params('scvelo') 
pk.plot_random_walks(
    seed=1121,
    n_sims=100,
    start_ixs={"bimod_oc_des": "ProgFP Development"},
    basis="umap",
    dpi=150,
)


# In[32]:


from cellrank.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata).compute_transition_matrix()


# In[33]:


combined_kernel = 0.8 * pk + 0.2 * ck
combined_kernel


# In[34]:


combined_kernel.compute_transition_matrix()


# In[35]:


scv.set_figure_params('scvelo') 
combined_kernel.plot_random_walks(
    seed=1121,
    n_sims=100,
    start_ixs={"bimod_oc_des": "ProgFP Development"},
    basis="umap",
    legend_loc="right",
    dpi=150,
)


# https://cellrank.readthedocs.io/en/stable/_images/100_cellrank_kernels.jpg

# # 3. Estimators

# In[36]:


from cellrank.estimators import GPCCA
g = GPCCA(pk)
print(g)


# ### Identify initial & terminal states

# In[37]:


g.fit(n_states=7, cluster_key="bimod_oc_des")
g.plot_macrostates(which="all")


# In[38]:


g.predict_terminal_states(method="top_n", n_states=6)
g.plot_macrostates(which="terminal")


# ### Compute fate probabilities and driver genes

# In[39]:


g.compute_fate_probabilities()
g.plot_fate_probabilities(legend_loc="right")


# In[40]:


g.plot_fate_probabilities(same_plot=False)


# In[41]:


cr.pl.circular_projection(adata, keys="orig.ident", legend_loc="right", title='', palette=oi_colors)


# In[42]:


mono_drivers = g.compute_lineage_drivers(lineages="NbM_2")
mono_drivers.head(10)


# In[43]:


mono_drivers = g.compute_lineage_drivers(lineages="DA Neuron projection")
mono_drivers.head(10)


# In[44]:


mono_drivers = g.compute_lineage_drivers(lineages="DA Synaptic assembly_1")
mono_drivers.head(10)


# In[45]:


mono_drivers = g.compute_lineage_drivers(lineages="DA Synaptic assembly_2")
mono_drivers.head(10)


# In[46]:


mono_drivers = g.compute_lineage_drivers(lineages="DA Synaptic modulation")
mono_drivers.head(10)


# In[47]:


df = mono_drivers.head(100)
df.to_csv('~/Desktop/XMAS_analysis/outputs/Dynamics/0729_cr_DA.SM_drivergenes_top100.csv', index=True)


# In[48]:


mono_drivers.head(100).index


# ### Visualize expression trends

# In[49]:


model = cr.models.GAMR(adata)


# In[50]:


cr.pl.gene_trends(
    adata,
    model=model,
    data_key="MAGIC_imputed_data",
    genes=["NTNG1", "GALNTL6", "CCSER1", "NKAIN2"],
    same_plot=True,
    ncols=2,
    time_key="palantir_pseudotime",
    hide_cells=True)


# In[51]:


cr.pl.gene_trends(
    adata,
    model=model,
    data_key="MAGIC_imputed_data",
    genes=["NTNG1", "GALNTL6", "CCSER1", "NKAIN2"],
    same_plot=True,
    ncols=2,
    time_key="palantir_pseudotime",
    hide_cells=True,
    legend_loc=None)


# In[52]:


cr.pl.heatmap(
    adata,
    model=model,
    data_key="MAGIC_imputed_data",
    genes=["NTNG1", "GALNTL6", "CCSER1", "NKAIN2"],
    lineages=["DA Neuron projection","DA Synaptic assembly_1", "DA Synaptic assembly_2", "DA Synaptic modulation"],
    time_key="palantir_pseudotime",
    cbar=False,
    show_all_genes=True,
)


# ### Aggregated fate probabilities

# In[53]:


states = ["ProgFP Development","NbM","DA Neuron projection","DA Synaptic assembly","DA Neurotransmitter release"]
sc.pl.embedding(
    adata, basis="umap", color="bimod_oc_des", groups=states, legend_loc="right"
)


# In[54]:


cr.pl.aggregate_fate_probabilities(
    adata,
    mode="violin",
    lineages=["DA Synaptic modulation"],
    cluster_key="bimod_oc_des",
    clusters=states,
)


# # 4. RNA velocity

# In[55]:


scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2


# In[56]:


import warnings

warnings.simplefilter("ignore", category=UserWarning)


# In[57]:


sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


# In[59]:


scv.tl.recover_dynamics(adata, n_jobs=5)
scv.tl.velocity(adata, mode="dynamical")


# In[60]:


scv.tl.velocity_graph(adata)


# In[61]:


adata


# In[62]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.velocity_embedding_grid(adata, basis='umap', color='bimod_oc_des')


# In[63]:


adata.write_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/s04_XMAS_velocyto_calculated.h5ad")


# In[ ]:


adata =  sc.read_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/s04_XMAS_velocyto_calculated.h5ad")


# ## 4.1 Velocity kernel

# In[64]:


vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()


# In[65]:


ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

combined_kernel_2 = 0.8 * vk + 0.2 * ck


# In[66]:


print(combined_kernel_2)


# In[67]:


scv.set_figure_params('scvelo') 
vk.plot_projection(basis="umap", color='orig.ident')


# In[68]:


scv.set_figure_params('scvelo') 
vk.plot_projection(basis="umap", color='palantir_pseudotime')


# In[69]:


vk.plot_random_walks(start_ixs={"bimod_oc_des": "ProgFP Development"}, max_iter=200, seed=1121)


# In[70]:


top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:10], ncols=5, frameon=False, color='orig.ident')


# In[71]:


scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters_BiMod_OC', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()


# In[72]:


scv.tl.rank_velocity_genes(adata, groupby='bimod_oc_des', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(10)


# In[259]:


# WRONG!!!
scv.tl.rank_velocity_genes(adata, groupby='bimod_oc_des', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(10)


# In[ ]:





# ## 4.2 Pseudotime

# In[73]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color='bimod_oc_des', add_outline=True, legend_loc='on data', size=120, title='',
                          legend_fontsize=11, legend_fontoutline=5,frameon=False )


# In[74]:


sc.pp.neighbors(adata, n_neighbors=100, n_pcs=40, use_rep='X_pca', method='gauss')
sc.tl.diffmap(adata)


# In[75]:


scv.tl.velocity_pseudotime(adata, vkey='velocity')


# In[76]:


Tn_mask = np.isin(adata.obs['bimod_oc_des'], ['ProgFP Development'])
min_stem_id = np.argmin(adata.obsm['X_diffmap'][Tn_mask, 1])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata.uns['iroot'] = root_id
scv.pl.scatter(adata,basis='umap',color=[root_id,'bimod_oc_des'],legend_loc='right')


# In[77]:


sc.set_figure_params(figsize=(5, 4))
sc.tl.dpt(adata)
adata.obsm['X_diffmap_'] = adata.obsm['X_diffmap'][:,1:]
sc.pl.embedding(adata,'diffmap_',color=['orig.ident','bimod_oc_des',])


# In[78]:


sc.tl.dpt(adata)
sc.pl.embedding(
    adata,
    basis="umap",
    color=["dpt_pseudotime","velocity_pseudotime", "palantir_pseudotime"],
    color_map="viridis_r",
)


# In[79]:


trajectory = ['ProgFP Development', 'NbM', 'DA Neuron projection', 
              'DA Synaptic assembly', 'DA Neurotransmitter release', 'DA Synaptic modulation']
trajectory_2 = ['DA Neuron projection', 'DA Synaptic assembly', 'DA Neurotransmitter release', 'DA Synaptic modulation']
mask = np.in1d(adata.obs["bimod_oc_des"], trajectory)
sc.pl.violin(
    adata[mask],
    keys=["dpt_pseudotime", "velocity_pseudotime", "palantir_pseudotime"],
    groupby="bimod_oc_des",
    rotation=-90,
    order=trajectory,
)
mask = np.in1d(adata.obs["bimod_oc_des"], trajectory_2)
sc.pl.violin(
    adata[mask],
    keys=["dpt_pseudotime","velocity_pseudotime", "palantir_pseudotime"],
    groupby="bimod_oc_des",
    rotation=-90,
    order=trajectory_2,
)


# ___

# ### 4.2.1 dpt_pseudotime

# In[80]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[83]:


scv.set_figure_params('scvelo') 
pk.plot_projection(color='dpt_pseudotime', cmap='viridis_r')


# In[84]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[85]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_BiMod_OC','bimod_oc_des'])


# In[86]:


pk.plot_projection(basis="umap", recompute=True, color=['bimod_oc_des'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[87]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='bimod_oc_des', n_convolve=100)
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='orig.ident', n_convolve=100)


# ### 4.2.2 velocity_pseudotime

# In[88]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="velocity_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[89]:


pk.plot_projection(color='velocity_pseudotime', cmap='viridis_r')


# In[90]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[91]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_BiMod_OC','bimod_oc_des'])


# In[92]:


pk.plot_projection(basis="umap", recompute=True, color=['bimod_oc_des'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[93]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='bimod_oc_des', n_convolve=100)
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='orig.ident', n_convolve=100)


# ### 4.2.3 palantir_pseudotime

# In[94]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="palantir_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[95]:


pk.plot_projection(color='palantir_pseudotime', cmap='viridis_r')


# In[96]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[97]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_BiMod_OC','bimod_oc_des'])


# In[98]:


pk.plot_projection(basis="umap", recompute=True, color=['bimod_oc_des'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[99]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='bimod_oc_des', n_convolve=100)
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='orig.ident', n_convolve=100)


# ### 4.2.4 latent time

# In[100]:


scv.tl.latent_time(adata, vkey='velocity', root_key='iroot', end_key='end')
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', basis='X_diffmap', size=80)


# In[101]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="latent_time")
pk.compute_transition_matrix()

print(pk)


# In[102]:


pk.plot_projection(color='latent_time', cmap='viridis_r')


# In[103]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[104]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_BiMod_OC','bimod_oc_des'])


# In[105]:


pk.plot_projection(basis="umap", recompute=True, color=['bimod_oc_des'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[106]:


pk.plot_projection(basis="umap", recompute=True, color=['bimod_oc_des'], legend_loc='right', palette='Paired', title='Velocity projection of DA lineage')


# In[107]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, color_map='rocket_r', sortby='latent_time', col_color='bimod_oc_des', n_convolve=200)
scv.pl.heatmap(adata, var_names=top_genes, color_map='rocket_r', sortby='latent_time', col_color='orig.ident', n_convolve=200)


# In[108]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, yticklabels=False, var_names=top_genes, color_map='rocket_r', 
               sortby='latent_time',figsize=(4,4), col_color='bimod_oc_des', n_convolve=300)
filtered_genes = [
    gene for gene in top_genes
    if not gene.startswith("RP") and
    "." not in gene and
    "-" not in gene and
    not gene.startswith("LINC")
]
scv.pl.heatmap(adata, yticklabels=False, var_names=filtered_genes, color_map='rocket_r', 
               sortby='latent_time',figsize=(4,4), col_color='bimod_oc_des', n_convolve=300)


# In[109]:


scv.pl.scatter(adata, color='latent_time', color_map='gnuplot_r', basis='umap', size=80,
              add_outline=True, legend_loc='on data',legend_fontsize=11, legend_fontoutline=2,frameon=False )


# In[110]:


scv.pl.scatter(adata, color='latent_time', color_map='gnuplot_r', basis='X_diffmap', size=80,
              add_outline=True, legend_loc='on data',legend_fontsize=11, legend_fontoutline=2,frameon=False )


# In[111]:


scv.tl.rank_velocity_genes(adata, groupby='bimod_oc_des', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(10)


# In[112]:


df.to_csv('/Users/gclyu07/Desktop/XMAS_analysis/outputs/Dynamics/0813_scvelo_drivers_bimod_oc_des.csv', index=False)


# In[ ]:





# In[113]:


adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='bimod_oc_des')


# In[114]:


scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[115]:


adata.write_h5ad("/Users/gclyu07/Desktop/XMAS_analysis/outputs/s04_4_2_4_0813_XMAS_velocyto_calculated.h5ad")


# _______

# In[116]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.scatter(adata, color=['orig.ident','bimod_oc_des',], color_map='gnuplot', basis='diffmap', size=70,
              add_outline=True, legend_loc='right',legend_fontsize=10, legend_fontoutline=2,frameon=False, title = ' ' )


# In[118]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot_r', basis='diffmap', size=80,
              add_outline=True, legend_loc='on data',legend_fontsize=11, legend_fontoutline=2,frameon=False, title = ' ' )


# In[119]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.scatter(adata, color=['bimod_oc',], basis='diffmap', size=70,
              add_outline=True, legend_loc='right',legend_fontsize=10, legend_fontoutline=2,frameon=False, title = ' ' )


# In[120]:


scv.pl.velocity_embedding_grid(adata, basis='umap', color='velocity_pseudotime')


# In[121]:


pk.plot_projection(basis="umap", recompute=True, color=['bimod_oc_des'], legend_loc='right', title=' ')


# In[ ]:




