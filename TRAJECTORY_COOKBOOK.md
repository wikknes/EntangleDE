# EntangleDE Trajectory Analysis Cookbook

This cookbook provides step-by-step guidance on using the quantum trajectory analysis functionality in EntangleDE. This feature uses quantum computing approaches to analyze single-cell RNA-seq data with time series information, providing enhanced trajectory inference capabilities.

## Table of Contents

1. [Introduction to Trajectory Analysis](#introduction-to-trajectory-analysis)
2. [Installation](#installation)
3. [Data Preparation](#data-preparation)
4. [Basic Usage](#basic-usage)
5. [Advanced Usage](#advanced-usage)
6. [Visualization](#visualization)
7. [Interpreting Results](#interpreting-results)
8. [Performance Considerations](#performance-considerations)
9. [Frequently Asked Questions](#frequently-asked-questions)

## Introduction to Trajectory Analysis

### What is Trajectory Analysis?

In single-cell RNA sequencing (scRNA-seq), trajectory analysis is the process of ordering cells along developmental or response pathways based on their gene expression profiles. Unlike static snapshots, trajectory analysis helps us understand:

- How cells transition between different states
- Developmental processes and cell differentiation
- Branching points where cells diverge into different cell types
- Gene expression dynamics along developmental pathways

### How Quantum Computing Enhances Trajectory Analysis

EntangleDE's trajectory analysis uses quantum computing in several key ways:

1. **Hamiltonian Embedding**: Cellular data is encoded as a quantum Hamiltonian, capturing the complex relationships and dynamics between cells.

2. **Quantum Clustering**: Uses quantum approximate optimization algorithms (QAOA) or quantum annealing to identify cell clusters more effectively.

3. **Quantum-Optimized Force Directed Graphs**: Constructs trajectory graphs with edge weights determined by quantum measurements.

These quantum approaches can reveal subtle patterns in the data that might be missed by classical methods, particularly for complex branching trajectories or noisy datasets.

## Installation

EntangleDE and its trajectory analysis functionality requires several Python packages. The core package depends on:

```bash
pip install numpy pandas matplotlib seaborn qiskit qiskit-aer scipy scikit-learn
```

For the trajectory analysis functionality, you'll also need:

```bash
pip install scanpy networkx
```

Optional: For quantum annealing capabilities (D-Wave backend):

```bash
pip install dwave-ocean-sdk
```

## Data Preparation

Trajectory analysis requires gene expression data and, optionally, time point information for each cell.

### Required Data Format

1. **Gene Expression Matrix**: Genes in rows, cells in columns
   - Normalized counts (e.g., log-normalized or TPM)
   - Dimensions: n_genes × n_cells

2. **Pseudotime/Time Series Data** (optional):
   - Vector of time/pseudotime values for each cell
   - Same length as the number of cells

### Example: Loading Data from Files

```python
import pandas as pd
import numpy as np
from src.pipeline import load_data

# Load expression data, pseudotime, and gene names
expression_data, pseudotime, gene_names = load_data(
    expression_file="path/to/expression_matrix.csv",
    pseudotime_file="path/to/pseudotime_values.csv",
    gene_names_file="path/to/gene_names.txt"
)
```

### Example: Using Synthetic Data

For testing or demonstration purposes, you can generate synthetic data with branching patterns:

```python
# Generate synthetic data with branching structure
from examples.trajectory_example import generate_branching_data

expr_data, pseudotime, branch_labels, gene_names = generate_branching_data(
    n_genes=100,   # Number of genes
    n_cells=200,   # Number of cells
    n_branches=3,  # Number of trajectory branches
    noise_level=0.1  # Amount of noise to add
)
```

## Basic Usage

The simplest way to run trajectory analysis is through the main pipeline:

```python
from src.pipeline import run_pipeline

results = run_pipeline(
    expression_data,  # Gene expression matrix (genes × cells)
    pseudotime,       # Pseudotime values for each cell
    gene_names,       # List of gene names
    n_components=20,  # Number of PCA components
    time_param=1.0,   # Time parameter for Hamiltonian
    n_measurements=1024,  # Number of quantum measurements
    output_dir="output",  # Directory to save results
    run_classical=True,   # Run classical analysis for comparison
    run_trajectory=True   # Enable trajectory analysis
)

# Access trajectory results
trajectory_results = results['trajectory']['results']
```

### Direct Trajectory Analysis

For more control, you can run trajectory analysis directly:

```python
from src.trajectory_analysis import quantum_trajectory_analysis

trajectory_results = quantum_trajectory_analysis(
    expression_data,  # Gene expression matrix (genes × cells)
    pseudotime,       # Pseudotime values for each cell
    gene_names,       # List of gene names
    n_components=20,  # Number of PCA components
    time_param=1.0,   # Time parameter for Hamiltonian
    n_measurements=1024,  # Number of quantum measurements
    quantum_backend='qiskit',  # 'qiskit' or 'dwave'
    n_clusters=5,     # Number of clusters for intermediate clustering
    output_dir="trajectory_output"  # Directory to save results
)
```

## Advanced Usage

For maximum flexibility, use the `QuantumTrajectoryAnalysis` class directly with AnnData objects:

```python
import scanpy as sc
from src.trajectory_analysis import QuantumTrajectoryAnalysis, load_as_anndata

# Convert data to AnnData format
adata = load_as_anndata(expression_data, pseudotime, gene_names)

# Create trajectory analyzer
analyzer = QuantumTrajectoryAnalysis(
    n_components=20,
    time_param=1.0,
    n_measurements=1024,
    quantum_backend='qiskit',  # or 'dwave' for quantum annealing
    n_neighbors=15  # Number of neighbors for graph construction
)

# Run trajectory analysis
results = analyzer.run_trajectory_analysis(
    adata,
    pseudotime,  # Optional: provide known pseudotime for comparison
    n_clusters=5  # Number of clusters for intermediate steps
)

# Compute quality metrics
metrics = analyzer.compute_trajectory_metrics(results['adata'], pseudotime)
print(f"Kendall's tau with true pseudotime: {metrics['kendall_tau']:.3f}")
print(f"Trajectory stability score: {metrics['stability']:.3f}")

# Plot results
analyzer.plot_trajectory(results['adata'], 
                        save_path="trajectory_plot.png")

# Plot gene dynamics along trajectory
genes_of_interest = ['Gene1', 'Gene2', 'Gene3']
analyzer.plot_gene_dynamics(results['adata'], 
                          genes_of_interest,
                          save_path="gene_dynamics.png")

# Export results
analyzer.export_trajectory("trajectory_results")
```

### Customizing the Quantum Backend

EntangleDE supports multiple quantum approaches:

#### Qiskit (Default)

Uses QAOA (Quantum Approximate Optimization Algorithm) for clustering cells:

```python
analyzer = QuantumTrajectoryAnalysis(
    quantum_backend='qiskit',
    n_measurements=1024  # More measurements = higher accuracy but slower
)
```

#### D-Wave Quantum Annealing

For larger problems, quantum annealing may be more efficient:

```python
analyzer = QuantumTrajectoryAnalysis(
    quantum_backend='dwave',
    # D-Wave specific settings are handled automatically
)
```

## Visualization

EntangleDE provides several visualization functions for trajectory analysis:

### 1. Trajectory Plot

Visualizes cells in a low-dimensional space, colored by clusters and pseudotime:

```python
analyzer.plot_trajectory(adata, 
                        embedding_key='X_umap',  # or 'X_pca'
                        save_path="trajectory.png")
```

### 2. Gene Expression Dynamics

Shows how gene expression changes along the trajectory:

```python
analyzer.plot_gene_dynamics(adata,
                          genes=['Gene1', 'Gene2', 'Gene3'],
                          save_path="gene_dynamics.png")
```

### 3. Custom Visualizations

You can create custom visualizations using the AnnData object:

```python
import matplotlib.pyplot as plt
import scanpy as sc

# Get results
adata = results['adata']
clusters = adata.obs['quantum_clusters']
pseudotime = adata.obs['quantum_pseudotime']

# Run UMAP if not already present
if 'X_umap' not in adata.obsm:
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)

# Create custom plot
plt.figure(figsize=(12, 5))

# Plot clusters
plt.subplot(1, 2, 1)
sc.pl.umap(adata, color='quantum_clusters', show=False)
plt.title('Cell Clusters')

# Plot pseudotime
plt.subplot(1, 2, 2)
sc.pl.umap(adata, color='quantum_pseudotime', show=False)
plt.title('Pseudotime')

plt.tight_layout()
plt.savefig("custom_trajectory_plot.png")
```

## Interpreting Results

The trajectory analysis produces several key results:

### 1. Quantum Pseudotime

A refined ordering of cells along their developmental trajectories:

```python
pseudotime = results['adata'].obs['quantum_pseudotime']
```

### 2. Cell Clusters

Identified cell groups that represent distinct cell states:

```python
clusters = results['adata'].obs['quantum_clusters']
```

### 3. Trajectory Graph

A force-directed graph representing cell relationships:

```python
graph = results['force_graph']
```

### 4. Quality Metrics

Metrics to evaluate trajectory inference quality:

```python
metrics = analyzer.compute_trajectory_metrics(results['adata'], true_pseudotime)
print(f"Kendall's tau: {metrics['kendall_tau']}")  # Correlation with true time
print(f"Stability: {metrics['stability']}")  # Robustness of trajectory
```

## Performance Considerations

### Dataset Size Recommendations

- **Small** (20-100 genes, 50-200 cells): Fast, good for testing
- **Medium** (100-500 genes, 200-1000 cells): Good balance of detail and performance
- **Large** (500+ genes, 1000+ cells): May require longer computation time

### Optimizing Performance

1. **Reduce dimensions**:
   ```python
   n_components=30  # Use fewer PCA components for large datasets
   ```

2. **Adjust measurement count**:
   ```python
   n_measurements=512  # Fewer measurements for faster runtime (default: 1024)
   ```

3. **For very large datasets**, pre-filter to highly variable genes:
   ```python
   # Use scanpy to select highly variable genes
   import scanpy as sc
   adata = load_as_anndata(expression_data)
   sc.pp.highly_variable_genes(adata, n_top_genes=500)
   adata = adata[:, adata.var.highly_variable]
   ```

## Frequently Asked Questions

### Q: When should I use quantum trajectory analysis vs. classical methods?

**A:** Quantum trajectory analysis is particularly beneficial for:
- Complex branching trajectories
- Noisy or sparse datasets
- When you need to identify subtle patterns
- When classical methods give inconsistent results

### Q: Do I need real quantum hardware to run this?

**A:** No, EntangleDE uses quantum simulators by default, which run on classical hardware. However, with the right API credentials, you can connect to real quantum computers (e.g., through IBM Quantum or D-Wave).

### Q: How do I interpret multiple branches in my trajectory?

**A:** Branches represent alternative cell fates or states. To analyze branches:
1. Look at the clustered cells in visualization
2. Examine branch-specific gene expression patterns
3. Correlate branches with known cell types or functions

### Q: My trajectory doesn't match my expected biology. What should I check?

**A:** Consider the following:
- Try adjusting `n_clusters` to better match expected cell states
- Ensure your data is properly normalized and scaled
- Check for batch effects in your data
- Verify if your genes capture the biological process of interest

### Q: Can I integrate this with existing Scanpy/Seurat workflows?

**A:** Yes! The `QuantumTrajectoryAnalysis` class works with AnnData objects, which are compatible with Scanpy. For Seurat (R), you would need to export/import data, but the workflow is compatible.

## Example Workflow

Here's a complete example showing a typical workflow:

```python
import numpy as np
import scanpy as sc
from src.trajectory_analysis import QuantumTrajectoryAnalysis, load_as_anndata

# 1. Load or generate data
expression_data = np.random.rand(100, 200)  # 100 genes, 200 cells
pseudotime = np.sort(np.random.rand(200))  # Random pseudotime
gene_names = [f"Gene_{i}" for i in range(100)]

# 2. Convert to AnnData
adata = load_as_anndata(expression_data, pseudotime, gene_names)

# 3. Initialize analyzer
analyzer = QuantumTrajectoryAnalysis(
    n_components=20,
    quantum_backend='qiskit'
)

# 4. Run analysis
results = analyzer.run_trajectory_analysis(adata)

# 5. Evaluate results
metrics = analyzer.compute_trajectory_metrics(results['adata'])
print(f"Stability score: {metrics['stability']}")

# 6. Visualize
analyzer.plot_trajectory(results['adata'], save_path="trajectory.png")

# 7. Find genes that change along trajectory
# Calculate variance along pseudotime
pt = results['adata'].obs['quantum_pseudotime'].values
var_genes = []
for i in range(expression_data.shape[0]):
    # Correlation with pseudotime
    corr = np.corrcoef(expression_data[i], pt)[0, 1]
    var_genes.append((gene_names[i], abs(corr)))

# Sort by correlation
var_genes.sort(key=lambda x: x[1], reverse=True)
print("Top variable genes along trajectory:")
for gene, corr in var_genes[:5]:
    print(f"{gene}: correlation = {corr:.3f}")

# 8. Plot top genes
top_genes = [g[0] for g in var_genes[:6]]
analyzer.plot_gene_dynamics(results['adata'], top_genes, save_path="top_genes.png")

# 9. Export results
analyzer.export_trajectory("trajectory_results")
```

This cookbook should help you get started with quantum trajectory analysis in EntangleDE. For more information, refer to the code documentation and API reference.