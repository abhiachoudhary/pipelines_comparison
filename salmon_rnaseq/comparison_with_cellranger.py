import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import time
def mprint(*variable, is_df=0): #my print
    for var in variable:
        if is_df == 1:
            print(var, '=\n', repr(eval(var)))
        else:
            print(var, '=', repr(eval(var)))

sample_type, sample_name = 'lung', 'GTEX-13N11-5030-SM-H5JDW.03.Lung.TST'
# sample_type, sample_name = 'heart', 'GTEX-1ICG6-5003-SM-GHS9A.01.Heart_Left_Ventricle.TST'

prefix_address = 'data/samples/'+sample_type+'/'+sample_name+'/'

#cellranger_org output
raw_expr_cr = sc.read_10x_h5(prefix_address+'cell_ranger_org/raw_gene_bc_matrices_h5.h5')
# raw_expr_cr.var_names_make_unique()
expr_cr = sc.read_10x_h5(prefix_address+'cell_ranger_org/filtered_gene_bc_matrices_h5.h5')
print(raw_expr_cr.shape) #(737280, 33694)
print(expr_cr.shape) #(6000, 33694)

#hubmap output
raw_expr_hm = sc.read_h5ad(prefix_address+'hubmap_result/raw_expr.h5ad')
expr_hm = sc.read_h5ad(prefix_address+'hubmap_result/expr.h5ad')
print(raw_expr_hm.shape) #(6107, 98000))
print(expr_hm.shape) #(6107, 60286))

#view before
cr_cell = expr_cr.obs.index.tolist()[:5]
cr_gene = expr_cr.var['gene_ids'].tolist()[:5]
hm_cell = expr_hm.obs.index.tolist()[:5]
hm_gene = expr_hm.var.index.tolist()[:5]
cr_gc_before = pd.DataFrame(list(zip(cr_cell, cr_gene)), index=['']*5, columns=['cell','gene'])
hm_gc_before = pd.DataFrame(list(zip(hm_cell, hm_gene)), index=['']*5, columns=['cell','gene'])
print(cr_gc_before)
print(hm_gc_before)

#mapping
expr_cr.obs.index = [x.split('-')[0] for x in expr_cr.obs.index.to_list()]
expr_hm.var.index = [x.split('.')[0] for x in expr_hm.var.index.to_list()]
expr_cr.var['hugo_symbol'] = expr_cr.var.index
expr_cr.var.index = expr_cr.var['gene_ids'].to_list()

#view after
cr_cell = expr_cr.obs.index.tolist()[:5]
cr_gene = expr_cr.var['gene_ids'].tolist()[:5]
hm_cell = expr_hm.obs.index.tolist()[:5]
hm_gene = expr_hm.var.index.tolist()[:5]
cr_gc_after = pd.DataFrame(list(zip(cr_cell, cr_gene)), index=['']*5, columns=['cell','gene'])
hm_gc_after = pd.DataFrame(list(zip(hm_cell, hm_gene)), index=['']*5, columns=['cell','gene'])
print(cr_gc_after)
print(hm_gc_after)

comp1 = expr_cr.copy()
comp2 = expr_hm.copy()
sample_names = ['cr', 'hm']

mprint('comp1')
mprint('comp2')

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata1 = comp1.copy()
adata2 = comp2.copy()

for adata in [adata1, adata2]:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], log1p=True, inplace=True)

# min_genes=200
# min_cells=3
# for adata in [adata1, adata2]:
#     sc.pp.filter_cells(adata, min_genes=min_genes)
#     sc.pp.filter_genes(adata, min_cells=min_cells)

# def more_filtering(adata):
#     adata = adata[adata.obs.n_genes_by_counts < 3000, :]
#     adata = adata[adata.obs.pct_counts_mt < 5, :]
#     return adata
    
# adata1 = more_filtering(adata1)
# adata2 = more_filtering(adata2)

for sample, adata in zip(sample_names, [adata1, adata2]):
    print(sample)
    adata.var_names_make_unique()
    ax = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts','pct_counts_mt'], 
                 jitter=0.4, multi_panel=True)

good_fig_size = (12,4)
nrows, ncols = 1, 2
fig = plt.figure(figsize=good_fig_size)
axes = [ fig.add_subplot(nrows, ncols, r * ncols + c + 1) for r in range(0, nrows) for c in range(0, ncols) ]
xvar, yvar = "log1p_total_counts", "log1p_n_genes_by_counts"
x_l, x_r, y_l, y_r = np.inf, -np.inf, np.inf, -np.inf
for adata in [adata1, adata2]:
    x_l = min(x_l, int(np.floor(adata.obs[xvar].min())))
    x_r = max(x_r, int(np.ceil(adata.obs[xvar].max())))
    y_l = min(y_l, int(np.floor(adata.obs[yvar].min())))
    y_r = max(y_r, int(np.ceil(adata.obs[yvar].max())))

for idx, (sample, adata) in enumerate(zip(sample_names,[adata1, adata2])):
    this_ax = axes[idx]
    sns.scatterplot(data=adata.obs, x=xvar,y=yvar, s=2, alpha=1, ax=this_ax, color='red')
    this_ax.set_xlim((x_l,x_r))
    this_ax.set_xlim((y_l,y_r))
    this_ax.set_title(sample)
    
    
#now we see the cells and genes overlap. Really good
mprint('len(np.unique(adata1.obs.index))')
mprint('len(np.unique(adata2.obs.index))')
common_cells = list(set(adata1.obs.index.to_list()) & set(adata2.obs.index.to_list()))
mprint('len(common_cells)')

mprint('len(np.unique(adata1.var.index))')
mprint('len(np.unique(adata2.var.index))')
common_genes = list(set(adata1.var.index.to_list()) & set(adata2.var.index.to_list()))
mprint('len(common_genes)')


df1 = adata1.to_df().loc[common_cells, common_genes]
df2 = adata2.to_df().loc[common_cells, common_genes]


#create a list with length of common_cells that contains the correlation between 
#gene expressions in expr_cr and expr_hm 
corr_vec_cw_prsn, corr_vec_cw_sprm = [], [] #cw = cellwise
for cell in common_cells:
    a1, a2 = df1.loc[cell,:], df2.loc[cell,:]
    a = pd.concat([a1, a2], axis=1)
    a1, a2 = a.iloc[:,0].to_list(), a.iloc[:,1].to_list()
    # corr_vec_cw_prsn.append(np.corrcoef(a1,a2)[0,1]) #pearson correlation
    corr_vec_cw_prsn.append(stats.pearsonr(a1,a2)[0]) #pearson correlation
    corr_vec_cw_sprm.append(stats.spearmanr(a1,a2)[0]) #spearman correlation

#create a list with length of common_genes that contains the correlation between 
#gene expressions in expr_cr and expr_hm 
corr_vec_gw_prsn, corr_vec_gw_sprm = [], [] #cw = cellwise
for gene in common_genes:
    a1, a2 = df1.loc[:,gene], df2.loc[:,gene]
    a = pd.concat([a1, a2], axis=1)
    a1, a2 = a.iloc[:,0].to_list(), a.iloc[:,1].to_list()
    # corr_vec_gw_prsn.append(np.corrcoef(a1,a2)[0,1]) #pearson correlation
    corr_vec_gw_prsn.append(stats.pearsonr(a1,a2)[0]) #pearson correlation
    corr_vec_gw_sprm.append(stats.spearmanr(a1,a2)[0]) #spearman correlation

plt.figure(figsize=(12, 5.5))
plt.subplot(221)
_ = plt.hist(corr_vec_cw_prsn, 50, density=0, facecolor='b', alpha=0.75, edgecolor='black')
plt.xlabel('Correlation (Pearson)')
plt.ylabel('Count (genes)')
plt.xlim(0,1)
plt.subplot(222)
_ = plt.hist(corr_vec_cw_sprm, 50, density=0, facecolor='b', alpha=0.75, edgecolor='black')
plt.xlabel('Correlation (Spearman)')
plt.ylabel('Count (genes)')
plt.xlim(0,1)

plt.subplot(223)
plt.text(0, 0, pd.Series(corr_vec_cw_prsn).describe().to_string(),)
plt.grid(False)
plt.axis('off')
plt.subplot(224)
plt.text(0, 0, pd.Series(corr_vec_cw_sprm).describe().to_string(),)
plt.grid(False)
plt.axis('off')

plt.figure(figsize=(12, 5.5))
plt.subplot(221)
_ = plt.hist(corr_vec_gw_prsn, 50, density=0, facecolor='b', alpha=0.75, edgecolor='black')
plt.xlabel('Correlation (Pearson)')
plt.ylabel('Count (cells)')
plt.xlim(0,1)
plt.subplot(222)
_ = plt.hist(corr_vec_gw_sprm, 50, density=0, facecolor='b', alpha=0.75, edgecolor='black')
plt.xlabel('Correlation (Spearman)')
plt.ylabel('Count (cells)')
plt.xlim(0,1)

plt.subplot(223)
plt.text(0, 0, pd.Series(corr_vec_gw_prsn).describe().to_string(),)
plt.grid(False)
plt.axis('off')
plt.subplot(224)
plt.text(0, 0, pd.Series(corr_vec_gw_sprm).describe().to_string(),)
plt.grid(False)
plt.axis('off')
