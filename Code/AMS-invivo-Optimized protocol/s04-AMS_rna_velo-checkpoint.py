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


# In[219]:


import os

# Path to datasets folder relative to your current directory
OS_path = "datasets"

# Full path to the .h5ad file
adata_path = os.path.join(OS_path, "05_neu_afterAnnotation.h5ad")

seed = 1121

# Optional: check if file exists first
print("File exists:", os.path.exists(adata_path))

# Read the h5ad file
adata = sc.read_h5ad(adata_path)

print(os.getcwd())

# Print object info to confirm load
print(adata)


# In[220]:


# Read the second h5ad file in the same directory
adata_u = sc.read_h5ad(os.path.join(OS_path, "s5_041725_neu_unspliced.h5ad"))

# Optional: check that it loaded correctly
print(adata_u)


# In[221]:


adata.X


# In[222]:


adata_u.X


# In[223]:


if np.all(adata_u.X == 0):
    print("All values in the matrix are zero.")
else:
    print("The matrix has non-zero values.")


# In[224]:


print(adata.var_names)
print(adata_u.var_names)


# In[225]:


overlapping_features = adata.var_names.intersection(adata_u.var_names)


# In[226]:


overlapping_features


# In[227]:


adata = adata[:, adata.var_names.isin(overlapping_features)]


# In[228]:


adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata_u.X


# In[229]:


adata


# In[230]:


adata.obs['neu.class'].value_counts()


# In[231]:


alist = {
    0: "Immature Neuron", 
    1: "DA", 
    2: "Glut",
    3: "GABA"
}
adata.obs['neu.class'] = (
adata.obs['neu.class']
.astype('category')
.map(alist)
)


# In[232]:


adata.obs['neu.class'].value_counts()


# In[233]:


alist = {
    0: "Immature Neuron",  # 0_NProg
    1: "DA_CALB1",
    2: "Glut",
    3: "GABA",
    4: "Glut",
    5: "Glut",
    6: "Glut",
    7: "DA_SOX6",
    8: "Immature Neuron",  # 8_NProg
    9: "Glut",
    10: "DA_SOX6",
    11: "Glut",
    12: "Glut",
    13: "DA_ONECUT2",
    14: "GABA",
    15: "Glut",
    16: "Glut",
    17: "Glut",
    18: "GABA"
}

adata.obs['neu.class_2'] = (
adata.obs['seurat_clusters_pca_RNA_2']
.astype('category')
.map(alist)
)


# In[234]:


adata.obs['neu.class_2'].value_counts()


# In[235]:


adata.obs['seurat_clusters_pca_RNA_2'] = adata.obs['seurat_clusters_pca_RNA_2'].astype('category')


# In[236]:


sc.pl.umap(adata, color = ['seurat_clusters_pca_RNA_2'])


# In[237]:


sc.pl.umap(adata, color = ['neu.class'])


# In[238]:


sc.pl.umap(adata, color = ['neu.class_2'])


# In[49]:


adata.write_h5ad(os.path.join(OS_path, "160725_RGS_velocyto_spliced_unspliced_corrected.h5ad"))


# In[241]:


adata_full = adata


# In[242]:


adata = adata_full[(adata_full.obs['neu.class'].isin(['DA','Immature Neuron']))& (adata_full.obs['seurat_clusters_pca_RNA_2'] != 8)]


# In[243]:


sc.set_figure_params(figsize=(4, 3))
sc.pl.umap(adata, color='orig.ident')
sc.pl.umap(adata, color = ['neu.class_2'])


# ___

# # 1. Palantir preprocessing

# In[244]:


sc.pp.neighbors(adata, n_neighbors=40, n_pcs=20, use_rep='X_pca', method='gauss')


# In[245]:


sc.tl.umap(adata)


# In[246]:


sc.set_figure_params(figsize=(5, 3))
sc.pl.umap(adata, color=['orig.ident','neu.class','neu.class_2'], add_outline=True, legend_loc='on data',
                          legend_fontsize=8, legend_fontoutline=2,frameon=False, title=['Original identities','Cell types'])


# In[247]:


dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)


# In[248]:


ms_data = palantir.utils.determine_multiscale_space(adata)


# In[249]:


imputed_X = palantir.utils.run_magic_imputation(adata)


# In[250]:


palantir.plot.plot_diffusion_components(adata)


# In[251]:


Tn_mask = np.isin(adata.obs['neu.class_2'], ['Immature Neuron'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask, 5])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata[[root_id]].obs_names


# In[252]:


palantir.plot.highlight_cells_on_umap(adata, pd.Series(["Immature Neuron"],
                index=["GTAATGAAGGTCATGC-3"]))


# In[253]:


Tn_mask = np.isin(adata.obs['neu.class_2'], ['DA_SOX6'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask,1])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata[[root_id]].obs_names


# In[254]:


Tn_mask = np.isin(adata.obs['neu.class_2'], ['DA_CALB1'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask,1])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata[[root_id]].obs_names


# In[255]:


Tn_mask = np.isin(adata.obs['neu.class_2'], ['DA_ONECUT2'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask,1])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata[[root_id]].obs_names


# In[256]:


start_cell = 'GTAATGAAGGTCATGC-3'
terminal_states = pd.Series(["DA_SOX6","DA_CALB1","DA_ONECUT2"],
                index=["ATCCCAAAGTTGCTCT-2", "GGAAATTGTAAGTACG-2", "GGCATCCGTTGACCTG-3"])


# In[257]:


palantir.plot.highlight_cells_on_umap(adata, terminal_states)


# In[258]:


pr_res = palantir.core.run_palantir(
    adata, start_cell, num_waypoints=500, terminal_states=terminal_states.index
)
pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]


# In[260]:


palantir.plot.plot_palantir_results(adata, s = 2.5)


# In[261]:


sc.pl.embedding(adata, basis="umap", color=["orig.ident","palantir_pseudotime","neu.class_2"], cmap="viridis_r")


# # 2. Kernels combination

# In[262]:


from cellrank.kernels import PseudotimeKernel
pk = PseudotimeKernel(adata, time_key="palantir_pseudotime")


# In[264]:


get_ipython().system('pip install ipywidgets')


# In[265]:


pk.compute_transition_matrix()


# In[266]:


scv.set_figure_params('scvelo') 
pk.plot_random_walks(
    seed=1121,
    n_sims=100,
    start_ixs={"neu.class_2": "Immature Neuron"},
    basis="umap",
    dpi=150,
)


# In[267]:


from cellrank.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata).compute_transition_matrix()


# In[268]:


combined_kernel = 0.8 * pk + 0.2 * ck
combined_kernel


# In[269]:


combined_kernel.compute_transition_matrix()


# In[270]:


scv.set_figure_params('scvelo') 
combined_kernel.plot_random_walks(
    seed=1121,
    n_sims=100,
    start_ixs={"neu.class_2": "Immature Neuron"},
    basis="umap",
    legend_loc="right",
    dpi=150,
)


# https://cellrank.readthedocs.io/en/stable/_images/100_cellrank_kernels.jpg

# # 3. Estimators

# In[271]:


from cellrank.estimators import GPCCA
g = GPCCA(pk)
print(g)


# ### Identify initial & terminal states

# In[272]:


g.fit(n_states=7, cluster_key="neu.class_2")
g.plot_macrostates(which="all")


# In[273]:


g.predict_terminal_states(method="top_n", n_states=6)
g.plot_macrostates(which="terminal")


# ### Compute fate probabilities and driver genes

# In[274]:


pip install petsc petsc4py


# In[275]:


g.compute_fate_probabilities()
g.plot_fate_probabilities(legend_loc="right")


# In[276]:


g.plot_fate_probabilities(same_plot=False)


# In[277]:


cr.pl.circular_projection(adata, keys="orig.ident", legend_loc="right", title='')


# In[278]:


mono_drivers = g.compute_lineage_drivers(lineages="DA_SOX6_1")
mono_drivers.head(10)


# In[279]:


mono_drivers = g.compute_lineage_drivers(lineages="DA_SOX6_2")
mono_drivers.head(10)


# In[280]:


mono_drivers = g.compute_lineage_drivers(lineages="DA_ONECUT2")
mono_drivers.head(10)

df = mono_drivers.head(100)
df.to_csv(os.path.join(OS_path, "0729_cr_DA.SM_drivergenes_top100.csv"), index=True)
# In[281]:


mono_drivers.head(100).index


# ### Visualize expression trends

# In[102]:


get_ipython().system('pip install rpy2')


# In[283]:


model = cr.models.GAMR(adata)


# In[284]:


cr.pl.gene_trends(
    adata,
    model=model,
    data_key="MAGIC_imputed_data",
    genes=["NTNG1", "GALNTL6", "CCSER1", "NKAIN2"],
    same_plot=True,
    ncols=2,
    time_key="palantir_pseudotime",
    hide_cells=True)


# In[285]:


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


# In[286]:


cr.pl.heatmap(
    adata,
    model=model,
    data_key="MAGIC_imputed_data",
    genes=["NTNG1", "GALNTL6", "CCSER1", "NKAIN2"],
    lineages=["DA_SOX6_1","DA_SOX6_2", "DA_CALB1_1", "DA_CALB1_2", "DA_ONECUT2"],
    time_key="palantir_pseudotime",
    cbar=False,
    show_all_genes=True,
)


# ### Aggregated fate probabilities

# In[287]:


states = ["Immature Neuron","DA_SOX6","DA_CALB1","DA_ONECUT2"]
sc.pl.embedding(
    adata, basis="umap", color="neu.class_2", groups=states, legend_loc="right"
)


# In[288]:


cr.pl.aggregate_fate_probabilities(
    adata,
    mode="violin",
    lineages=["DA_SOX6_1"],
    cluster_key="neu.class_2",
    clusters=states,
)


# In[289]:


cr.pl.aggregate_fate_probabilities(
    adata,
    mode="violin",
    lineages=["DA_ONECUT2"],
    cluster_key="neu.class_2",
    clusters=states,
)


# # 4. RNA velocity

# In[293]:


scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2


# In[294]:


import warnings

warnings.simplefilter("ignore", category=UserWarning)


# In[295]:


import multiprocessing
print(multiprocessing.cpu_count())


# In[296]:


sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


# In[297]:


# Preprocessing (if not yet done)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# Fit dynamics
scv.tl.recover_dynamics(adata, n_jobs=20)

# Compute velocity
scv.tl.velocity(adata, mode="dynamical")

# Optional: compute velocity graph and visualize
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')


# In[298]:


adata


# In[299]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.velocity_embedding_grid(adata, basis='umap', color='neu.class_2')


# In[208]:


adata.write_h5ad(os.path.join(OS_path, "160725_RGS_velocyto_calculated.h5ad"))

adata = sc.read_h5ad(os.path.join(OS_path, "160725_RGS_velocyto_calculated.h5ad"))
# ## 4.1 Velocity kernel

# In[300]:


vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()


# In[301]:


ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

combined_kernel_2 = 0.8 * vk + 0.2 * ck


# In[302]:


print(combined_kernel_2)


# In[303]:


scv.set_figure_params('scvelo') 
vk.plot_projection(basis="umap", color='orig.ident')


# In[304]:


scv.set_figure_params('scvelo') 
vk.plot_projection(basis="umap", color='palantir_pseudotime')


# In[305]:


vk.plot_random_walks(start_ixs={"neu.class_2": "Immature Neuron"}, max_iter=200, seed=1121)


# In[327]:


import matplotlib.pyplot as plt


# In[319]:


# Compute mean expression per gene (columns = genes)
# If sparse matrix, convert to dense or use appropriate method
if hasattr(adata.X, "toarray"):
    mean_expression = np.array(adata.X.toarray()).mean(axis=0)
else:
    mean_expression = np.array(adata.X).mean(axis=0)

# Create a Series with gene names as index
mean_exp_series = pd.Series(mean_expression, index=adata.var_names)

# Sort descending and get top 10
top10_expressed = mean_exp_series.sort_values(ascending=False).head(10).index.tolist()

print("Top 10 most expressed genes:", top10_expressed)


# In[321]:


print("Top 10 most dynamic genes:",top_genes[:10].tolist())


# In[328]:


top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:10], ncols=5, color='neu.class_2')


# In[322]:


scv.pl.scatter(adata, basis=top10_expressed, ncols=5, color='neu.class_2')


# In[315]:


genes_of_interest = ["TH", "NURR1", "EN1", "SOX6", "CALB1", "ONECUT2", "PITX3", "OTX2", "SLC6A3", "SLC18A2", "LMX1A", "DDC", "FOXA2", "KCNJ6"]

# Filter genes that are actually in the dataset
genes_to_plot = [gene for gene in genes_of_interest if gene in adata.var_names]

# Plot phase portraits
scv.pl.scatter(adata, genes_to_plot, color='neu.class_2', ncols = 3)


# In[133]:


scv.tl.rank_velocity_genes(adata, groupby='neu.class_2', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()


# In[135]:


scv.tl.rank_velocity_genes(adata, groupby='neu.class_2', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(10)


# In[ ]:





# ## 4.2 Pseudotime

# In[136]:


sc.set_figure_params(figsize=(8, 6))
sc.pl.umap(adata, color='neu.class_2', add_outline=True, legend_loc='on data', size=120, title='',
                          legend_fontsize=11, legend_fontoutline=5,frameon=False )


# In[137]:


sc.pp.neighbors(adata, n_neighbors=100, n_pcs=40, use_rep='X_pca', method='gauss')
sc.tl.diffmap(adata)


# In[138]:


scv.tl.velocity_pseudotime(adata, vkey='velocity')


# In[139]:


Tn_mask = np.isin(adata.obs['neu.class_2'], ['Immature Neuron'])
min_stem_id = np.argmin(adata.obsm['X_pca'][Tn_mask, 5])
root_id = np.arange(len(Tn_mask))[Tn_mask][min_stem_id]
adata.uns['iroot'] = root_id
scv.pl.scatter(adata,basis='umap',color=[root_id,'neu.class_2'],legend_loc='right')


# In[140]:


sc.set_figure_params(figsize=(5, 4))
sc.tl.dpt(adata)
adata.obsm['X_diffmap_'] = adata.obsm['X_diffmap'][:,1:]
sc.pl.embedding(adata,'diffmap_',color=['orig.ident','neu.class_2',])


# In[141]:


sc.tl.dpt(adata)
sc.pl.embedding(
    adata,
    basis="umap",
    color=["dpt_pseudotime","velocity_pseudotime", "palantir_pseudotime"],
    color_map="viridis_r",
)


# In[142]:


trajectory = ['Immature Neuron',"DA_CALB1","DA_ONECUT2","DA_SOX6"]
trajectory_2 = ["DA_CALB1","DA_ONECUT2","DA_SOX6"]
mask = np.in1d(adata.obs["neu.class_2"], trajectory)
sc.pl.violin(
    adata[mask],
    keys=["dpt_pseudotime", "velocity_pseudotime", "palantir_pseudotime"],
    groupby="neu.class_2",
    rotation=-90,
    order=trajectory,
)
mask = np.in1d(adata.obs["neu.class_2"], trajectory_2)
sc.pl.violin(
    adata[mask],
    keys=["dpt_pseudotime","velocity_pseudotime", "palantir_pseudotime"],
    groupby="neu.class_2",
    rotation=-90,
    order=trajectory_2,
)


# ___

# ### 4.2.1 dpt_pseudotime

# In[143]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[144]:


scv.set_figure_params('scvelo') 
pk.plot_projection(color='dpt_pseudotime', cmap='viridis_r')


# In[145]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[146]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_pca_RNA_2','neu.class_2'])


# In[147]:


pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[148]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='neu.class_2', n_convolve=100)


# In[149]:


pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[150]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='neu.class_2', n_convolve=100)


# ### 4.2.2 velocity_pseudotime

# In[151]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="velocity_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[152]:


pk.plot_projection(color='velocity_pseudotime', cmap='viridis_r')


# In[153]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[154]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_pca_RNA_2','neu.class_2'])


# In[155]:


pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[156]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='neu.class_2', n_convolve=100)


# ### 4.2.3 palantir_pseudotime

# In[157]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="palantir_pseudotime")
pk.compute_transition_matrix()

print(pk)


# In[158]:


pk.plot_projection(color='palantir_pseudotime', cmap='viridis_r')


# In[159]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='orig.ident')


# In[160]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_pca_RNA_2','neu.class_2'])


# In[161]:


pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], legend_loc='right', palette='Spectral', title='Velocity projection of DA lineage')


# In[162]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='dpt_pseudotime', col_color='neu.class_2', n_convolve=100)


# ### 4.2.4 latent time

# In[163]:


scv.tl.latent_time(adata, vkey='velocity', root_key='iroot', end_key='end')
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', basis='X_diffmap', size=80)


# In[164]:


pk = cr.kernels.PseudotimeKernel(adata, time_key="latent_time")
pk.compute_transition_matrix()

print(pk)


# In[165]:


pk.plot_projection(color='latent_time', cmap='viridis_r')


# In[166]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="diffmap", recompute=True, color='neu.class_2')


# In[167]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','seurat_clusters_pca_RNA_2','neu.class_2'])


# In[168]:


pk.plot_projection(basis="umap", recompute=True, color=['orig.ident','neu.class_2','latent_time'], 
                   legend_loc='right', palette='Spectral', title=['Original identity','Velocity projection of DA lineage', 'Latent time'])


# In[169]:


import seaborn
seaborn.color_palette("Paired")


# In[170]:


scv.set_figure_params('scvelo') 
pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], 
                   legend_loc='right', palette='Paired', title=[''])


# In[171]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, color_map='rocket_r', sortby='latent_time', col_color='neu.class_2',palette='Paired', n_convolve=200)


# In[173]:


pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], legend_loc='right', palette='Paired', title='Velocity projection of DA lineage')


# In[174]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, color_map='rocket_r', sortby='latent_time', col_color='neu.class_2', n_convolve=200)


# In[175]:


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, yticklabels=False, var_names=top_genes, color_map='rocket_r', 
               sortby='latent_time',figsize=(6,4), col_color='neu.class_2', n_convolve=300)
filtered_genes = [
    gene for gene in top_genes
    if not gene.startswith("RP") and
    "." not in gene and
    "-" not in gene and
    not gene.startswith("LINC")
]
scv.pl.heatmap(adata, yticklabels=False, var_names=filtered_genes, color_map='rocket_r', 
               sortby='latent_time',figsize=(4,4), col_color='neu.class_2', n_convolve=200)


# In[176]:


scv.pl.scatter(adata, color='latent_time', color_map='rocket_r', basis='umap', size=80,
              add_outline=True, legend_loc='on data',legend_fontsize=11, legend_fontoutline=2,frameon=False )


# In[177]:


scv.pl.scatter(adata, color='latent_time', color_map='gnuplot_r', basis='X_diffmap', size=80,
              add_outline=True, legend_loc='on data',legend_fontsize=11, legend_fontoutline=2,frameon=False )


# In[178]:


scv.tl.rank_velocity_genes(adata, groupby='neu.class_2', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(10)


# In[179]:


scv.tl.rank_velocity_genes(adata, groupby='neu.class_2', min_corr=.3)

df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])

# Filter function
def filter_genes(gene):
    return not gene.startswith("ENSG") and "." not in gene and "-" not in gene and not gene.startswith("LINC")

# Apply the filter and remove non-matching genes
filtered_df = df.applymap(lambda gene: gene if filter_genes(gene) else None)

filtered_df = filtered_df.apply(lambda col: col.dropna().reset_index(drop=True))

filtered_df.head(10)


# In[180]:


df.to_csv(os.path.join(OS_path, "1607RGS_scvelo_drivers_bimod_oc_des.csv"), index=False)


# In[ ]:





# In[182]:


adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='neu.class_2')


# In[ ]:


scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[245]:


adata.write_h5ad("/Users/gclyu07/Desktop/R_data/outputs/050425_f_XMAS_velocyto_calculated.h5ad")


# _______

# In[246]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.scatter(adata, color=['orig.ident','neu.class_2',], color_map='gnuplot', basis='diffmap', size=70,
              add_outline=True, legend_loc='right',legend_fontsize=10, legend_fontoutline=2,frameon=False, title = ' ' )


# In[247]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot_r', basis='diffmap', size=80,
              add_outline=True, legend_loc='on data',legend_fontsize=11, legend_fontoutline=2,frameon=False, title = ' ' )


# In[249]:


scv.settings.verbosity = 0  # only show errors, no hints/information
scv.set_figure_params('scvelo')  # for beautified visualization
scv.pl.scatter(adata, color=['neu.class_2',], basis='diffmap', size=70,
              add_outline=True, legend_loc='right',legend_fontsize=10, legend_fontoutline=2,frameon=False, title = ' ' )


# In[250]:


scv.pl.velocity_embedding_grid(adata, basis='umap', color='velocity_pseudotime')


# In[251]:


pk.plot_projection(basis="umap", recompute=True, color=['neu.class_2'], legend_loc='right', title=' ')


# In[ ]:




