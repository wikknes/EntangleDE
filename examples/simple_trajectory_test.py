import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.trajectory_analysis import quantum_trajectory_analysis

# Generate a simple trajectory dataset
n_genes = 20
n_cells = 50

# Create pseudotime from 0 to 1
pseudotime = np.linspace(0, 1, n_cells)

# Initialize expression data with some noise
expression_data = np.random.rand(n_genes, n_cells) * 0.5

# Add time-dependent patterns to first 5 genes
for i in range(5):
    if i % 2 == 0:  # Linear pattern
        expression_data[i, :] = 3 * pseudotime + np.random.randn(n_cells) * 0.1
    else:  # Quadratic pattern
        expression_data[i, :] = 3 * (pseudotime ** 2) + np.random.randn(n_cells) * 0.1

# Create gene names
gene_names = [f"Gene_{i}" for i in range(n_genes)]

# Create output directory
os.makedirs("simple_test_results", exist_ok=True)

# Run trajectory analysis
print("Running trajectory analysis...")
trajectory_results = quantum_trajectory_analysis(
    expression_data,
    pseudotime,
    gene_names,
    n_components=10,
    time_param=1.0,
    n_measurements=1024,
    quantum_backend='qiskit',
    output_dir="simple_test_results"
)

print("Analysis complete! Results saved to simple_test_results/")