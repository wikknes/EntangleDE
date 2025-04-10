import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.pipeline import load_data, run_pipeline

# Create a simple synthetic dataset
def generate_synthetic_data(n_genes=100, n_cells=200, n_diff_genes=10):
    # Create pseudotime
    pseudotime = np.linspace(0, 1, n_cells)
    
    # Initialize expression data
    expression_data = np.random.rand(n_genes, n_cells) * 0.5
    
    # Create different gene expression patterns along pseudotime
    for i in range(n_diff_genes):
        if i % 3 == 0:  # Linear increasing
            expression_data[i, :] = i + pseudotime * 5 + np.random.randn(n_cells) * 0.2
        elif i % 3 == 1:  # Sigmoidal
            midpoint = 0.5
            steepness = 10
            expression_data[i, :] = 5 / (1 + np.exp(-steepness * (pseudotime - midpoint))) + np.random.randn(n_cells) * 0.2
        else:  # Peak in the middle
            expression_data[i, :] = 3 * np.exp(-((pseudotime - 0.5) ** 2) / 0.05) + np.random.randn(n_cells) * 0.2
    
    # Create gene names
    gene_names = []
    for i in range(n_genes):
        if i < n_diff_genes:
            gene_names.append(f"DIFF_GENE_{i}")
        else:
            gene_names.append(f"STABLE_GENE_{i}")
    
    return expression_data, pseudotime, gene_names

# Generate data with a smaller dataset for testing
expr_data, pseudotime, gene_names = generate_synthetic_data(n_genes=20, n_cells=50, n_diff_genes=5)

# Create output directory
os.makedirs("example_results", exist_ok=True)

# Save synthetic data to files
pd.DataFrame(expr_data, index=gene_names).to_csv("example_results/synthetic_expression.csv")
pd.DataFrame({"pseudotime": pseudotime}).to_csv("example_results/synthetic_pseudotime.csv", index=False)
with open("example_results/synthetic_gene_names.txt", "w") as f:
    for gene in gene_names:
        f.write(f"{gene}\n")

# Plot some example gene patterns
plt.figure(figsize=(12, 8))
for i in range(9):
    plt.subplot(3, 3, i+1)
    plt.plot(pseudotime, expr_data[i, :])
    plt.title(gene_names[i])
    plt.xlabel("Pseudotime")
    plt.ylabel("Expression")
plt.tight_layout()
plt.savefig("example_results/example_patterns.png")
plt.close()

# Run the analysis
print("\nRunning Quantum Differential Expression Analysis...")
results = run_pipeline(
    expr_data, 
    pseudotime,
    gene_names,
    n_components=20,
    time_param=1.0,
    n_measurements=1024,
    output_dir="example_results"
)

# Print the top identified genes
top_indices = results['quantum']['diff_results']['diff_genes_indices'][:10]
top_genes = [gene_names[i] for i in top_indices]
print("\nTop differentially expressed genes identified:")
for i, gene in enumerate(top_genes):
    print(f"{i+1}. {gene}")

# Check how many true differentially expressed genes were identified
true_diff_genes = [g for g in top_genes if "DIFF" in g]
print(f"\nAccuracy: {len(true_diff_genes)}/10 true differential genes in the top 10")

# Calculate overall accuracy percentage for all differentially expressed genes
top_20_indices = results['quantum']['diff_results']['diff_genes_indices'][:20]
top_20_genes = [gene_names[i] for i in top_20_indices]
true_diff_genes_20 = [g for g in top_20_genes if "DIFF" in g]
accuracy = len(true_diff_genes_20) / 10 * 100
print(f"Overall accuracy: {accuracy}% of true differential genes identified in top 20")
print(f"\nResults have been saved to the 'example_results' directory")
