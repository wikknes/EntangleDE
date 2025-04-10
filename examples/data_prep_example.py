import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

"""
This script demonstrates how to prepare real single-cell RNA-seq data 
for use with the QDTA pipeline. It shows how to:
1. Load a real scRNA-seq dataset
2. Preprocess the data
3. Calculate pseudotime
4. Export the data in the correct format for QDTA
"""

# Create output directory
os.makedirs("example_data", exist_ok=True)

# Option 1: Use a built-in example dataset from scanpy
print("Loading example data...")
adata = sc.datasets.paul15()
print(f"Loaded data with {adata.shape[0]} cells and {adata.shape[1]} genes")

# Basic preprocessing
print("Preprocessing data...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
adata = adata[:, adata.var.highly_variable]
print(f"After preprocessing: {adata.shape[0]} cells and {adata.shape[1]} genes")

# Perform PCA and KNN graph construction
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

# UMAP for visualization
sc.tl.umap(adata)

# Calculate pseudotime using diffusion pseudotime
print("Calculating pseudotime...")
# For this example, we'll use diffusion pseudotime
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_dcs=10, n_branchings=0)  # Not considering branching for simplicity

# Plot UMAP with pseudotime coloring
plt.figure(figsize=(10, 8))
sc.pl.umap(adata, color=['dpt_pseudotime'], cmap='viridis', show=False)
plt.savefig("example_data/pseudotime_umap.png")
plt.close()

# Export data for QDTA
print("Exporting data for QDTA...")

# Extract the expression matrix (genes x cells)
expr_matrix = adata.X.T  # Transpose to get genes x cells

# Convert to DataFrame with gene names
gene_names = adata.var_names.tolist()
expr_df = pd.DataFrame(expr_matrix, index=gene_names)

# Export expression matrix
expr_df.to_csv("example_data/paul15_expression.csv")

# Export pseudotime
pseudotime_df = pd.DataFrame({
    'cell': adata.obs_names,
    'pseudotime': adata.obs['dpt_pseudotime']
})
pseudotime_df.to_csv("example_data/paul15_pseudotime.csv", index=False)

# Export gene names
with open("example_data/paul15_gene_names.txt", "w") as f:
    for gene in gene_names:
        f.write(f"{gene}\n")

print("Data preparation complete. Files saved to 'example_data' directory.")
print("You can now run QDTA with these files using:")
print("python main.py \\
  --expression example_data/paul15_expression.csv \\
  --pseudotime example_data/paul15_pseudotime.csv \\
  --genes example_data/paul15_gene_names.txt \\
  --output paul15_results")
