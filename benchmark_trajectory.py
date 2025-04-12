#!/usr/bin/env python
'''python benchmark_trajectory.py --size small
  python benchmark_trajectory.py --size medium
  python benchmark_trajectory.py --size large'''
import numpy as np
import pandas as pd
import argparse
import time
import os
import psutil
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, silhouette_score
from scipy.stats import spearmanr, kendalltau

from src.trajectory_analysis import QuantumTrajectoryAnalysis, load_as_anndata

def generate_synthetic_trajectory_data(n_genes=100, n_cells=200, n_branches=3, noise_level=0.1):
    """
    Create a synthetic dataset with a branching trajectory structure.
    
    Args:
        n_genes (int): Number of genes
        n_cells (int): Number of cells
        n_branches (int): Number of trajectory branches
        noise_level (float): Level of noise to add
        
    Returns:
        tuple: (expression_data, true_time, branch_labels, gene_names)
    """
    print(f"Generating synthetic trajectory data with {n_genes} genes, {n_cells} cells, {n_branches} branches...")
    
    # Create pseudotime from 0 to 1
    pseudotime = np.random.rand(n_cells)
    
    # Sort for convenience
    pseudotime = np.sort(pseudotime)
    
    # Create branch assignments - start with branch 0, 
    # then randomly assign to branches after a branching point
    branch_point = 0.3  # Branching occurs at pseudotime 0.3
    early_cells = pseudotime < branch_point
    late_cells = ~early_cells
    
    branch_labels = np.zeros(n_cells, dtype=int)
    # Assign random branches to late cells
    branch_labels[late_cells] = np.random.randint(0, n_branches, size=np.sum(late_cells))
    
    # Initialize expression data
    expression_data = np.zeros((n_genes, n_cells))
    
    # Create different gene expression patterns
    gene_index = 0
    
    # 1. Early genes (before branch point, then decrease)
    n_early_genes = n_genes // 4
    for i in range(n_early_genes):
        # Expression increases until branch point, then decreases
        pattern = np.zeros(n_cells)
        for j, t in enumerate(pseudotime):
            if t < branch_point:
                pattern[j] = t / branch_point
            else:
                pattern[j] = 1.0 - (t - branch_point) / (1.0 - branch_point)
        
        # Add noise
        pattern += np.random.normal(0, noise_level, n_cells)
        pattern = np.clip(pattern, 0, 1)
        
        expression_data[gene_index] = pattern
        gene_index += 1
    
    # 2. Branch-specific genes
    n_branch_genes = n_genes // 4
    for branch in range(n_branches):
        for i in range(n_branch_genes // n_branches):
            # Expression increases only in the specific branch
            pattern = np.zeros(n_cells)
            for j, t in enumerate(pseudotime):
                if t >= branch_point and branch_labels[j] == branch:
                    pattern[j] = (t - branch_point) / (1.0 - branch_point)
            
            # Add noise
            pattern += np.random.normal(0, noise_level, n_cells)
            pattern = np.clip(pattern, 0, 1)
            
            expression_data[gene_index] = pattern
            gene_index += 1
    
    # 3. Late genes (increase after branch point regardless of branch)
    n_late_genes = n_genes // 4
    for i in range(n_late_genes):
        # Expression increases after branch point
        pattern = np.zeros(n_cells)
        for j, t in enumerate(pseudotime):
            if t >= branch_point:
                pattern[j] = (t - branch_point) / (1.0 - branch_point)
        
        # Add noise
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
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    
    return expression_data, pseudotime, branch_labels, gene_names

def run_quantum_trajectory(expression_data, pseudotime, gene_names, config):
    """
    Run quantum trajectory analysis with the given configuration.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray): True pseudotime values
        gene_names (list): List of gene names
        config (dict): Configuration parameters
        
    Returns:
        dict: Analysis results and metrics
    """
    print(f"Running quantum trajectory analysis with backend: {config['quantum_backend']}")
    
    # Convert to AnnData
    adata = load_as_anndata(expression_data, pseudotime, gene_names)
    
    # Create analyzer
    analyzer = QuantumTrajectoryAnalysis(
        n_components=config['n_components'],
        time_param=config['time_param'],
        n_measurements=config['n_measurements'],
        quantum_backend=config['quantum_backend'],
        n_neighbors=config['n_neighbors']
    )
    
    # Start timing
    start_time = time.time()
    
    # Run analysis
    results = analyzer.run_trajectory_analysis(
        adata, 
        pseudotime, 
        n_clusters=config['n_clusters']
    )
    
    # End timing
    end_time = time.time()
    execution_time = end_time - start_time
    
    return {
        'results': results,
        'execution_time': execution_time,
        'analyzer': analyzer
    }

def run_scanpy_trajectory(expression_data, pseudotime, gene_names):
    """
    Run Scanpy's classical trajectory analysis for comparison.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray): True pseudotime values
        gene_names (list): List of gene names
        
    Returns:
        dict: Analysis results and metrics
    """
    print("Running classical Scanpy trajectory analysis for comparison")
    
    # Convert to AnnData
    adata = load_as_anndata(expression_data, pseudotime, gene_names)
    
    # Start timing
    start_time = time.time()
    
    # Basic preprocessing
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(1000, adata.shape[1]))
    adata = adata[:, adata.var.highly_variable]
    
    # PCA and neighborhood graph
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    
    # Run UMAP for visualization
    sc.tl.umap(adata)
    
    # Run Leiden clustering
    sc.tl.leiden(adata)
    
    # Run PAGA for trajectory inference
    sc.tl.paga(adata, groups='leiden')
    
    # Run diffusion pseudotime (DPT)
    sc.tl.diffmap(adata)
    
    # Set the root cell to be the one with lowest true pseudotime
    root_cell_idx = np.argmin(pseudotime)
    adata.uns['iroot'] = root_cell_idx
    
    # Calculate DPT
    sc.tl.dpt(adata)
    
    # End timing
    end_time = time.time()
    execution_time = end_time - start_time
    
    return {
        'adata': adata,
        'execution_time': execution_time
    }

def evaluate_trajectory(quantum_results, scanpy_results, true_pseudotime, branch_labels):
    """
    Evaluate trajectory inference results against ground truth.
    
    Args:
        quantum_results (dict): Quantum trajectory analysis results
        scanpy_results (dict): Scanpy trajectory analysis results
        true_pseudotime (np.ndarray): True pseudotime values
        branch_labels (np.ndarray): True branch assignments
        
    Returns:
        dict: Evaluation metrics
    """
    print("Evaluating trajectory inference results")
    
    metrics = {}
    
    # Get quantum results
    q_adata = quantum_results['results']['adata']
    q_pseudotime = quantum_results['results']['quantum_pseudotime']
    q_clusters = quantum_results['results']['quantum_clusters']
    
    # Get Scanpy results
    s_adata = scanpy_results['adata']
    s_pseudotime = s_adata.obs['dpt_pseudotime'].values
    s_clusters = s_adata.obs['leiden'].astype(int).values
    
    # 1. Pseudotime correlation
    q_tau, q_p = kendalltau(q_pseudotime, true_pseudotime)
    s_tau, s_p = kendalltau(s_pseudotime, true_pseudotime)
    
    metrics['pseudotime_correlation'] = {
        'quantum_kendall_tau': q_tau,
        'quantum_p_value': q_p,
        'scanpy_kendall_tau': s_tau,
        'scanpy_p_value': s_p
    }
    
    # 2. Clustering accuracy (compared to branches)
    q_ari = adjusted_rand_score(branch_labels, q_clusters)
    s_ari = adjusted_rand_score(branch_labels, s_clusters)
    
    metrics['clustering_accuracy'] = {
        'quantum_ari': q_ari,
        'scanpy_ari': s_ari
    }
    
    # 3. Silhouette score (internal validation)
    try:
        q_silhouette = silhouette_score(q_adata.obsm['X_pca'], q_clusters)
    except:
        q_silhouette = 0
    
    try:
        s_silhouette = silhouette_score(s_adata.obsm['X_pca'], s_clusters)
    except:
        s_silhouette = 0
    
    metrics['silhouette_score'] = {
        'quantum': q_silhouette,
        'scanpy': s_silhouette
    }
    
    # 4. Execution time comparison
    metrics['execution_time'] = {
        'quantum': quantum_results['execution_time'],
        'scanpy': scanpy_results['execution_time']
    }
    
    return metrics

def visualize_results(quantum_results, scanpy_results, metrics, true_pseudotime, 
                     branch_labels, output_dir):
    """
    Visualize trajectory analysis results.
    
    Args:
        quantum_results (dict): Quantum trajectory analysis results
        scanpy_results (dict): Scanpy trajectory analysis results
        metrics (dict): Evaluation metrics
        true_pseudotime (np.ndarray): True pseudotime values
        branch_labels (np.ndarray): True branch assignments
        output_dir (str): Directory to save visualizations
    """
    print("Visualizing results")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get data
    q_adata = quantum_results['results']['adata']
    q_pseudotime = quantum_results['results']['quantum_pseudotime']
    q_clusters = quantum_results['results']['quantum_clusters']
    
    s_adata = scanpy_results['adata']
    s_pseudotime = s_adata.obs['dpt_pseudotime'].values
    s_clusters = s_adata.obs['leiden'].astype(int).values
    
    # 1. Pseudotime comparison
    plt.figure(figsize=(15, 5))
    
    plt.subplot(1, 3, 1)
    plt.scatter(true_pseudotime, q_pseudotime, alpha=0.7)
    plt.xlabel('True Pseudotime')
    plt.ylabel('Quantum Pseudotime')
    plt.title(f"Quantum (τ={metrics['pseudotime_correlation']['quantum_kendall_tau']:.3f})")
    
    plt.subplot(1, 3, 2)
    plt.scatter(true_pseudotime, s_pseudotime, alpha=0.7)
    plt.xlabel('True Pseudotime')
    plt.ylabel('Scanpy Pseudotime')
    plt.title(f"Scanpy (τ={metrics['pseudotime_correlation']['scanpy_kendall_tau']:.3f})")
    
    plt.subplot(1, 3, 3)
    plt.scatter(q_pseudotime, s_pseudotime, alpha=0.7)
    plt.xlabel('Quantum Pseudotime')
    plt.ylabel('Scanpy Pseudotime')
    corr, _ = spearmanr(q_pseudotime, s_pseudotime)
    plt.title(f"Q vs S (ρ={corr:.3f})")
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pseudotime_comparison.png", dpi=300)
    plt.close()
    
    # 2. UMAP visualizations
    # If UMAP not already calculated for quantum, use the scanpy UMAP
    if 'X_umap' not in q_adata.obsm:
        q_adata.obsm['X_umap'] = s_adata.obsm['X_umap'].copy()
    
    umap_coords = q_adata.obsm['X_umap']
    
    # True branches
    plt.figure(figsize=(15, 12))
    
    plt.subplot(2, 2, 1)
    scatter = plt.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                         c=branch_labels, cmap='tab10', s=50, alpha=0.7)
    plt.colorbar(scatter, label='True Branches')
    plt.title('True Trajectory Branches')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    
    # True pseudotime
    plt.subplot(2, 2, 2)
    scatter = plt.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                         c=true_pseudotime, cmap='viridis', s=50, alpha=0.7)
    plt.colorbar(scatter, label='Pseudotime')
    plt.title('True Pseudotime')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    
    # Quantum clusters
    plt.subplot(2, 2, 3)
    scatter = plt.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                         c=q_clusters, cmap='tab10', s=50, alpha=0.7)
    plt.colorbar(scatter, label='Clusters')
    plt.title(f'Quantum Clusters (ARI={metrics["clustering_accuracy"]["quantum_ari"]:.3f})')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    
    # Quantum pseudotime
    plt.subplot(2, 2, 4)
    scatter = plt.scatter(umap_coords[:, 0], umap_coords[:, 1], 
                         c=q_pseudotime, cmap='viridis', s=50, alpha=0.7)
    plt.colorbar(scatter, label='Pseudotime')
    plt.title(f'Quantum Pseudotime (τ={metrics["pseudotime_correlation"]["quantum_kendall_tau"]:.3f})')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/trajectory_visualization.png", dpi=300)
    plt.close()
    
    # 3. Performance comparison
    plt.figure(figsize=(12, 6))
    
    methods = ['Quantum', 'Scanpy']
    
    # Execution time
    plt.subplot(1, 2, 1)
    times = [metrics['execution_time']['quantum'], metrics['execution_time']['scanpy']]
    plt.bar(methods, times)
    plt.ylabel('Execution Time (seconds)')
    plt.title('Performance Comparison')
    for i, v in enumerate(times):
        plt.text(i, v + 0.1, f"{v:.2f}s", ha='center')
    
    # Pseudotime accuracy
    plt.subplot(1, 2, 2)
    accuracy = [metrics['pseudotime_correlation']['quantum_kendall_tau'], 
               metrics['pseudotime_correlation']['scanpy_kendall_tau']]
    plt.bar(methods, accuracy)
    plt.ylabel("Kendall's τ with true pseudotime")
    plt.title('Pseudotime Accuracy')
    for i, v in enumerate(accuracy):
        plt.text(i, max(0, v) + 0.01, f"{v:.3f}", ha='center')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/performance_comparison.png", dpi=300)
    plt.close()
    
    # 4. Generate graph visualization if available
    try:
        analyzer = quantum_results['analyzer']
        force_graph = quantum_results['results']['force_graph']
        
        # Plot force-directed graph
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(force_graph)
        
        # Color nodes by pseudotime
        node_colors = [q_pseudotime[i] for i in force_graph.nodes()]
        
        # Draw the graph
        nx.draw_networkx_nodes(force_graph, pos, node_color=node_colors, 
                              cmap='viridis', node_size=50, alpha=0.8)
        nx.draw_networkx_edges(force_graph, pos, alpha=0.2)
        
        plt.title('Quantum Force-Directed Trajectory Graph')
        plt.axis('off')
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap='viridis')
        sm.set_array(node_colors)
        plt.colorbar(sm, label='Pseudotime')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/trajectory_graph.png", dpi=300)
        plt.close()
    except Exception as e:
        print(f"Error generating graph visualization: {e}")

def save_metrics(metrics, output_dir, size):
    """
    Save metrics to CSV file.
    
    Args:
        metrics (dict): Evaluation metrics
        output_dir (str): Directory to save metrics
        size (str): Dataset size
    """
    # Create flattened metrics for CSV
    flat_metrics = {
        'dataset_size': size,
        'quantum_execution_time': metrics['execution_time']['quantum'],
        'scanpy_execution_time': metrics['execution_time']['scanpy'],
        'quantum_kendall_tau': metrics['pseudotime_correlation']['quantum_kendall_tau'],
        'scanpy_kendall_tau': metrics['pseudotime_correlation']['scanpy_kendall_tau'],
        'quantum_ari': metrics['clustering_accuracy']['quantum_ari'],
        'scanpy_ari': metrics['clustering_accuracy']['scanpy_ari'],
        'quantum_silhouette': metrics['silhouette_score']['quantum'],
        'scanpy_silhouette': metrics['silhouette_score']['scanpy']
    }
    
    # Convert to DataFrame
    df = pd.DataFrame([flat_metrics])
    
    # Save to CSV
    df.to_csv(f"{output_dir}/{size}_trajectory_metrics.csv", index=False)

def run_benchmark(size='small', output_dir='benchmark_trajectory_results'):
    """
    Run benchmark for trajectory analysis.
    
    Args:
        size (str): Dataset size ('small', 'medium', 'large')
        output_dir (str): Directory to save results
        
    Returns:
        dict: Benchmark results
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Configure dataset based on size
    if size == 'small':
        n_genes, n_cells, n_branches = 20, 50, 2
        config = {
            'n_components': 10,
            'time_param': 1.0,
            'n_measurements': 1024,
            'quantum_backend': 'qiskit',
            'n_neighbors': 10,
            'n_clusters': 3
        }
    elif size == 'medium':
        n_genes, n_cells, n_branches = 100, 200, 3
        config = {
            'n_components': 20,
            'time_param': 1.0,
            'n_measurements': 1024,
            'quantum_backend': 'qiskit',
            'n_neighbors': 15,
            'n_clusters': 5
        }
    elif size == 'large':
        n_genes, n_cells, n_branches = 500, 500, 4
        config = {
            'n_components': 30,
            'time_param': 1.0,
            'n_measurements': 2048,
            'quantum_backend': 'qiskit',
            'n_neighbors': 20,
            'n_clusters': 6
        }
    else:
        raise ValueError(f"Invalid size: {size}")
    
    # Generate synthetic data
    expr_data, pseudotime, branch_labels, gene_names = generate_synthetic_trajectory_data(
        n_genes, n_cells, n_branches, noise_level=0.1
    )
    
    # Save synthetic data
    size_output_dir = os.path.join(output_dir, size)
    os.makedirs(size_output_dir, exist_ok=True)
    
    pd.DataFrame(expr_data, index=gene_names).to_csv(f"{size_output_dir}/synthetic_expression.csv")
    pd.DataFrame({"pseudotime": pseudotime}).to_csv(f"{size_output_dir}/synthetic_pseudotime.csv", index=False)
    pd.DataFrame({"branch": branch_labels}).to_csv(f"{size_output_dir}/synthetic_branches.csv", index=False)
    
    # Monitor memory
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss / (1024 * 1024)  # MB
    
    # Run quantum trajectory analysis
    quantum_results = run_quantum_trajectory(expr_data, pseudotime, gene_names, config)
    
    # Run Scanpy trajectory analysis for comparison
    scanpy_results = run_scanpy_trajectory(expr_data, pseudotime, gene_names)
    
    # Calculate peak memory usage
    peak_memory = process.memory_info().rss / (1024 * 1024)  # MB
    memory_used = peak_memory - initial_memory
    
    # Evaluate results
    metrics = evaluate_trajectory(quantum_results, scanpy_results, pseudotime, branch_labels)
    metrics['memory_usage'] = memory_used
    
    # Visualize results
    visualize_results(quantum_results, scanpy_results, metrics, pseudotime, 
                     branch_labels, size_output_dir)
    
    # Save metrics
    save_metrics(metrics, output_dir, size)
    
    # Print results
    print("\n===== TRAJECTORY BENCHMARK RESULTS =====")
    print(f"Dataset: {size} ({n_genes} genes, {n_cells} cells, {n_branches} branches)")
    print(f"Quantum execution time: {metrics['execution_time']['quantum']:.2f}s")
    print(f"Scanpy execution time: {metrics['execution_time']['scanpy']:.2f}s")
    print(f"Memory usage: {memory_used:.1f}MB")
    
    print("\nPseudotime accuracy:")
    print(f"  Quantum: Kendall's τ = {metrics['pseudotime_correlation']['quantum_kendall_tau']:.3f}")
    print(f"  Scanpy: Kendall's τ = {metrics['pseudotime_correlation']['scanpy_kendall_tau']:.3f}")
    
    print("\nClustering accuracy:")
    print(f"  Quantum: ARI = {metrics['clustering_accuracy']['quantum_ari']:.3f}")
    print(f"  Scanpy: ARI = {metrics['clustering_accuracy']['scanpy_ari']:.3f}")
    
    return {
        'quantum_results': quantum_results,
        'scanpy_results': scanpy_results,
        'metrics': metrics
    }

def main():
    parser = argparse.ArgumentParser(description='Run QDTA trajectory analysis benchmarking')
    parser.add_argument('--size', choices=['small', 'medium', 'large'], default='small',
                        help='Dataset size to benchmark (default: small)')
    parser.add_argument('--output', default='benchmark_trajectory_results',
                        help='Output directory for benchmark results (default: benchmark_trajectory_results)')
    args = parser.parse_args()
    
    run_benchmark(args.size, args.output)

if __name__ == "__main__":
    main()