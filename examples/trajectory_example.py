import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import sys
import scanpy as sc
from sklearn.preprocessing import StandardScaler

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.trajectory_analysis import quantum_trajectory_analysis, QuantumTrajectoryAnalysis
from src.pipeline import load_data, run_pipeline

# Create a synthetic branching trajectory dataset
def generate_branching_data(n_genes=100, n_cells=200, n_branches=3, noise_level=0.1):
    """Generate synthetic data with branching trajectory."""
    
    # Create pseudotime from 0 to 1
    pseudotime = np.random.rand(n_cells)
    pseudotime = np.sort(pseudotime)
    
    # Create branch assignments
    branch_point = 0.3  # Branching occurs at pseudotime 0.3
    early_cells = pseudotime < branch_point
    late_cells = ~early_cells
    
    branch_labels = np.zeros(n_cells, dtype=int)
    branch_labels[late_cells] = np.random.randint(0, n_branches, size=np.sum(late_cells))
    
    # Initialize expression data
    expression_data = np.zeros((n_genes, n_cells))
    
    # Create different gene expression patterns
    gene_index = 0
    
    # 1. Early genes (before branch point, then decrease)
    n_early_genes = n_genes // 4
    for i in range(n_early_genes):
        pattern = np.zeros(n_cells)
        for j, t in enumerate(pseudotime):
            if t < branch_point:
                pattern[j] = t / branch_point
            else:
                pattern[j] = 1.0 - (t - branch_point) / (1.0 - branch_point)
        
        pattern += np.random.normal(0, noise_level, n_cells)
        pattern = np.clip(pattern, 0, 1)
        expression_data[gene_index] = pattern
        gene_index += 1
    
    # 2. Branch-specific genes
    n_branch_genes = n_genes // 4
    for branch in range(n_branches):
        for i in range(n_branch_genes // n_branches):
            pattern = np.zeros(n_cells)
            for j, t in enumerate(pseudotime):
                if t >= branch_point and branch_labels[j] == branch:
                    pattern[j] = (t - branch_point) / (1.0 - branch_point)
            
            pattern += np.random.normal(0, noise_level, n_cells)
            pattern = np.clip(pattern, 0, 1)
            expression_data[gene_index] = pattern
            gene_index += 1
    
    # 3. Late genes (increase after branch point regardless of branch)
    n_late_genes = n_genes // 4
    for i in range(n_late_genes):
        pattern = np.zeros(n_cells)
        for j, t in enumerate(pseudotime):
            if t >= branch_point:
                pattern[j] = (t - branch_point) / (1.0 - branch_point)
        
        pattern += np.random.normal(0, noise_level, n_cells)
        pattern = np.clip(pattern, 0, 1)
        expression_data[gene_index] = pattern
        gene_index += 1
    
    # 4. Fill remaining genes with random noise
    while gene_index < n_genes:
        expression_data[gene_index] = np.random.normal(0.5, noise_level, n_cells)
        expression_data[gene_index] = np.clip(expression_data[gene_index], 0, 1)
        gene_index += 1
    
    # Create gene names
    gene_names = []
    gene_index = 0
    
    # Early genes
    for i in range(n_early_genes):
        gene_names.append(f"Early_Gene_{i}")
        gene_index += 1
    
    # Branch-specific genes
    for branch in range(n_branches):
        for i in range(n_branch_genes // n_branches):
            gene_names.append(f"Branch{branch}_Gene_{i}")
            gene_index += 1
    
    # Late genes
    for i in range(n_late_genes):
        gene_names.append(f"Late_Gene_{i}")
        gene_index += 1
    
    # Random genes
    while gene_index < n_genes:
        gene_names.append(f"Random_Gene_{gene_index}")
        gene_index += 1
    
    return expression_data, pseudotime, branch_labels, gene_names

# Generate synthetic data with branching structure
print("Generating synthetic branching data...")
expr_data, pseudotime, branch_labels, gene_names = generate_branching_data(
    n_genes=50, n_cells=100, n_branches=3
)

# Create output directory
os.makedirs("example_trajectory_results", exist_ok=True)

# Save synthetic data
pd.DataFrame(expr_data, index=gene_names).to_csv("example_trajectory_results/trajectory_expression.csv")
pd.DataFrame({"pseudotime": pseudotime}).to_csv("example_trajectory_results/trajectory_pseudotime.csv", index=False)
pd.DataFrame({"branch": branch_labels}).to_csv("example_trajectory_results/trajectory_branches.csv", index=False)

# Plot some example gene expression patterns
plt.figure(figsize=(15, 10))

# 1. Plot branch structure
plt.subplot(2, 3, 1)
# Create a 2D representation for visualization
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(StandardScaler().fit_transform(expr_data.T))

plt.scatter(pca_result[:, 0], pca_result[:, 1], c=branch_labels, cmap='tab10', s=50, alpha=0.7)
plt.title('Branching Structure (PCA)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.colorbar(label='Branch')

# 2. Plot pseudotime
plt.subplot(2, 3, 2)
plt.scatter(pca_result[:, 0], pca_result[:, 1], c=pseudotime, cmap='viridis', s=50, alpha=0.7)
plt.title('Pseudotime')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.colorbar(label='Pseudotime')

# 3. Plot examples of early genes
plt.subplot(2, 3, 3)
for i in range(3):
    idx = i
    plt.plot(pseudotime, expr_data[idx], label=gene_names[idx])
plt.title('Early Genes')
plt.xlabel('Pseudotime')
plt.ylabel('Expression')
plt.legend()

# 4. Plot examples of branch-specific genes
plt.subplot(2, 3, 4)
for i in range(3):
    idx = expr_data.shape[0] // 4 + i
    plt.plot(pseudotime, expr_data[idx], label=gene_names[idx])
plt.title('Branch 0 Genes')
plt.xlabel('Pseudotime')
plt.ylabel('Expression')
plt.legend()

# 5. Plot examples of late genes
plt.subplot(2, 3, 5)
for i in range(3):
    idx = expr_data.shape[0] // 2 + i
    plt.plot(pseudotime, expr_data[idx], label=gene_names[idx])
plt.title('Late Genes')
plt.xlabel('Pseudotime')
plt.ylabel('Expression')
plt.legend()

# 6. Plot examples of random genes
plt.subplot(2, 3, 6)
for i in range(3):
    idx = expr_data.shape[0] - 4 + i
    plt.plot(pseudotime, expr_data[idx], label=gene_names[idx])
plt.title('Random Genes')
plt.xlabel('Pseudotime')
plt.ylabel('Expression')
plt.legend()

plt.tight_layout()
plt.savefig("example_trajectory_results/gene_patterns.png")
plt.close()

# Method 1: Run through the main pipeline with trajectory analysis
print("\n==== Running Trajectory Analysis through Pipeline ====")
results = run_pipeline(
    expr_data, 
    pseudotime,
    gene_names,
    n_components=20,
    time_param=1.0,
    n_measurements=1024,
    output_dir="example_trajectory_results/pipeline_output",
    run_classical=True,
    run_trajectory=True  # Enable trajectory analysis
)

# Method 2: Run trajectory analysis directly with more options
print("\n==== Running Direct Trajectory Analysis ====")
trajectory_results = quantum_trajectory_analysis(
    expr_data,
    pseudotime,
    gene_names,
    n_components=20,
    time_param=1.0,
    n_measurements=1024,
    quantum_backend='qiskit',
    n_clusters=4,
    output_dir="example_trajectory_results/direct_output"
)

# Method 3: Use the QuantumTrajectoryAnalysis class directly with AnnData
print("\n==== Running Advanced Trajectory Analysis with AnnData ====")

# Convert to AnnData
from src.trajectory_analysis import load_as_anndata
adata = load_as_anndata(expr_data, pseudotime, gene_names)

# Create analyzer
analyzer = QuantumTrajectoryAnalysis(
    n_components=20,
    time_param=1.0,
    n_measurements=1024,
    quantum_backend='qiskit',
    n_neighbors=15
)

# Run analysis
advanced_results = analyzer.run_trajectory_analysis(adata, pseudotime, n_clusters=4)

# Compute metrics against ground truth
metrics = analyzer.compute_trajectory_metrics(advanced_results['adata'], pseudotime)
print(f"Kendall's tau with true pseudotime: {metrics['kendall_tau']:.3f}")
print(f"Trajectory stability score: {metrics['stability']:.3f}")

# Plot results
fig = analyzer.plot_trajectory(advanced_results['adata'], 
                             save_path="example_trajectory_results/advanced_trajectory.png")

# Plot expression of top branch-specific genes
branch_genes = [gene for gene in gene_names if 'Branch' in gene][:6]
gene_fig = analyzer.plot_gene_dynamics(advanced_results['adata'], branch_genes,
                                     save_path="example_trajectory_results/advanced_gene_dynamics.png")

# Export results to file
analyzer.export_trajectory("example_trajectory_results/advanced_output")

# Summarize execution time
print("\n==== Performance Summary ====")
if 'trajectory' in results:
    pipeline_time = results['trajectory']['execution_time']
    print(f"Pipeline trajectory analysis: {pipeline_time:.2f} seconds")

direct_time = trajectory_results['evaluation_metrics'].get('execution_time', 0)
print(f"Direct trajectory analysis: {direct_time:.2f} seconds")

print("\nResults have been saved to the 'example_trajectory_results' directory")