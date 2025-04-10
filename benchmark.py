#!/usr/bin/env python
'''python benchmark.py --size small
  python benchmark.py --size medium
  python benchmark.py --size large'''
import numpy as np
import pandas as pd
import argparse
import time
import os
import psutil
import matplotlib.pyplot as plt
from src.pipeline import load_data, run_pipeline

def generate_synthetic_data(n_genes=100, n_cells=200, n_diff_genes=10):
    """
    Create a synthetic dataset with genes showing different expression patterns.
    """
    # Create pseudotime
    pseudotime = np.linspace(0, 1, n_cells)
    
    # Initialize expression data
    expression_data = np.random.rand(n_genes, n_cells) * 0.5
    
    # Create different gene expression patterns along pseudotime
    for i in range(n_diff_genes):
        pattern_type = i % 3
        
        if pattern_type == 0:  # Linear increasing
            expression_data[i, :] = i + pseudotime * 5 + np.random.randn(n_cells) * 0.2
        elif pattern_type == 1:  # Sigmoidal
            midpoint = 0.5
            steepness = 10
            expression_data[i, :] = 5 / (1 + np.exp(-steepness * (pseudotime - midpoint))) + np.random.randn(n_cells) * 0.2
        else:  # Bell-curve
            expression_data[i, :] = 3 * np.exp(-((pseudotime - 0.5) ** 2) / 0.05) + np.random.randn(n_cells) * 0.2
    
    # Create gene names
    gene_names = []
    for i in range(n_genes):
        if i < n_diff_genes:
            gene_names.append(f"DIFF_GENE_{i}")
        else:
            gene_names.append(f"STABLE_GENE_{i}")
    
    return expression_data, pseudotime, gene_names

def calculate_metrics(results, gene_names, diff_count):
    """Calculate performance metrics."""
    # Get top differentially expressed genes
    quantum_indices = results['quantum']['diff_results']['diff_genes_indices']
    classical_indices = results['classical']['diff_genes_indices']
    
    # Calculate top-5, top-10, top-20 accuracy
    metrics = {}
    
    for n in [5, 10, 20]:
        n = min(n, diff_count)  # Don't exceed number of diff genes
        
        # Quantum accuracy
        q_top_n = quantum_indices[:n]
        q_correct = sum(1 for idx in q_top_n if "DIFF" in gene_names[idx])
        q_accuracy = q_correct / min(n, diff_count) * 100
        
        # Classical accuracy
        c_top_n = classical_indices[:n]
        c_correct = sum(1 for idx in c_top_n if "DIFF" in gene_names[idx])
        c_accuracy = c_correct / min(n, diff_count) * 100
        
        # Overlap
        overlap = len(set(q_top_n).intersection(set(c_top_n)))
        overlap_pct = overlap / n * 100
        
        metrics[f'top_{n}'] = {
            'quantum_accuracy': q_accuracy,
            'classical_accuracy': c_accuracy,
            'overlap': overlap_pct
        }
    
    return metrics

def run_benchmark(size='small', output_dir='benchmark_results'):
    """Run benchmark for specified dataset size."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Configure dataset based on size
    if size == 'small':
        n_genes, n_cells, n_diff_genes = 20, 50, 5
    elif size == 'medium':
        n_genes, n_cells, n_diff_genes = 100, 200, 10
    elif size == 'large':
        n_genes, n_cells, n_diff_genes = 500, 500, 30
    else:
        raise ValueError(f"Invalid size: {size}")
    
    # Generate synthetic data
    print(f"Generating {size} dataset ({n_genes} genes, {n_cells} cells, {n_diff_genes} diff genes)...")
    expr_data, pseudotime, gene_names = generate_synthetic_data(n_genes, n_cells, n_diff_genes)
    
    # Configure parameters
    n_components = min(20, n_genes // 2)
    
    # Monitor memory
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss / (1024 * 1024)  # MB
    
    # Run pipeline
    print(f"Running pipeline on {size} dataset...")
    output_subdir = os.path.join(output_dir, size)
    start_time = time.time()
    results = run_pipeline(
        expr_data, 
        pseudotime,
        gene_names,
        n_components=n_components,
        time_param=1.0,
        n_measurements=1024,
        output_dir=output_subdir
    )
    end_time = time.time()
    
    # Calculate peak memory usage
    peak_memory = process.memory_info().rss / (1024 * 1024)  # MB
    memory_used = peak_memory - initial_memory
    
    # Calculate metrics
    metrics = calculate_metrics(results, gene_names, n_diff_genes)
    
    # Print results
    print("\n===== BENCHMARK RESULTS =====")
    print(f"Dataset: {size} ({n_genes} genes, {n_cells} cells, {n_diff_genes} diff genes)")
    print(f"Quantum runtime: {results['quantum']['execution_time']:.2f}s")
    print(f"Classical runtime: {results['classical']['execution_time']:.2f}s")
    print(f"Memory usage: {memory_used:.1f}MB")
    
    print("\nAccuracy:")
    for n in [5, 10, 20]:
        if n > n_diff_genes:
            continue
        print(f"  Top-{n}: Quantum={metrics[f'top_{n}']['quantum_accuracy']:.1f}%, " + 
              f"Classical={metrics[f'top_{n}']['classical_accuracy']:.1f}%, " +
              f"Overlap={metrics[f'top_{n}']['overlap']:.1f}%")
    
    # Save metrics to CSV
    results_df = pd.DataFrame({
        'Metric': [
            'Quantum Runtime (s)',
            'Classical Runtime (s)',
            'Memory Usage (MB)',
            'Top-5 Quantum Accuracy (%)',
            'Top-5 Classical Accuracy (%)',
            'Top-5 Overlap (%)',
            'Top-10 Quantum Accuracy (%)',
            'Top-10 Classical Accuracy (%)',
            'Top-10 Overlap (%)',
            'Top-20 Quantum Accuracy (%)',
            'Top-20 Classical Accuracy (%)',
            'Top-20 Overlap (%)',
        ],
        'Value': [
            results['quantum']['execution_time'],
            results['classical']['execution_time'],
            memory_used,
            metrics.get('top_5', {}).get('quantum_accuracy', 'N/A'),
            metrics.get('top_5', {}).get('classical_accuracy', 'N/A'),
            metrics.get('top_5', {}).get('overlap', 'N/A'),
            metrics.get('top_10', {}).get('quantum_accuracy', 'N/A'),
            metrics.get('top_10', {}).get('classical_accuracy', 'N/A'),
            metrics.get('top_10', {}).get('overlap', 'N/A'),
            metrics.get('top_20', {}).get('quantum_accuracy', 'N/A'),
            metrics.get('top_20', {}).get('classical_accuracy', 'N/A'),
            metrics.get('top_20', {}).get('overlap', 'N/A'),
        ]
    })
    
    results_df.to_csv(os.path.join(output_dir, f"{size}_benchmark_results.csv"), index=False)
    
    return results, metrics

def main():
    parser = argparse.ArgumentParser(description='Run QDTA benchmarking')
    parser.add_argument('--size', choices=['small', 'medium', 'large'], default='small',
                        help='Dataset size to benchmark (default: small)')
    parser.add_argument('--output', default='benchmark_results',
                        help='Output directory for benchmark results (default: benchmark_results)')
    args = parser.parse_args()
    
    run_benchmark(args.size, args.output)

if __name__ == "__main__":
    main()