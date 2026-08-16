#!/usr/bin/env python3
"""Build a small synthetic scRNA-seq count matrix to demo/testing scEssentials."""
import numpy as np
import anndata as ad
import scipy.sparse as sp

rng = np.random.default_rng(42)

# 5000 genes, 300 cells
n_genes, n_cells = 5000, 300

# housekeeping-like stable high genes (will be detectable & stable)
stable = ['ACTB', 'GAPDH', 'EEF1A1', 'GAPDHS', 'B2M', 'TUBA1B', 'TUBB', 'RPL13A',
          'RPLP0', 'RPS18', 'PPIA', 'HPRT1', 'COX7C', 'ATP5F1A', 'PGK1', 'OAZ1',
          'RPL27', 'RPS13', 'HNRNPA1', 'EEF2', 'LDHA', 'ENO1', 'TPI1', 'ALDOA',
          'YWHAZ', 'SDHA', 'UBC', 'CYC1', 'NDUFA1', 'GPI']
other = [f'GENE{i:04d}' for i in range(n_genes - len(stable))]
genes = stable + other

# base log-normal expression per gene
mu = rng.normal(0.3, 0.6, n_genes)
mu = np.clip(mu, 0.0, 1.2)          # keep non-essential genes modestly expressed
mu[:len(stable)] = 2.5               # essential genes high
# add mild over-dispersion / noise to everything but keep essentials stable
expr_var = np.ones(n_genes)
expr_var[:len(stable)] = 0.2         # essential genes are LOW-variability

counts = np.zeros((n_genes, n_cells))
for g in range(n_genes):
    disp = expr_var[g]
    lam = np.exp(mu[g]) * rng.gamma(1.0 / disp, disp, n_cells)
    counts[g] = rng.poisson(lam)
counts = counts.astype(np.float32)

X = sp.csr_matrix(counts.T)  # cells x genes
adata = ad.AnnData(X)
adata.var_names = np.array(genes, dtype='U32')
adata.obs_names = [f'cell{i:03d}' for i in range(n_cells)]
adata.write_h5ad('data_demo.h5ad')

print('wrote data_demo.h5ad', counts.shape)
print('stable genes present:', len(set(stable) & set(adata.var_names)))
