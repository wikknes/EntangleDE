import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os

from src.quantum_gene_analysis import quantum_differential_analysis
from src.classical_benchmark import classical_differential_analysis, compare_methods
from src.trajectory_analysis import quantum_trajectory_analysis, QuantumTrajectoryAnalysis

def load_data(expression_file, pseudotime_file=None, gene_names_file=None):
    """
    Load gene expression and pseudotime data.
    
    Args:
        expression_file (str): Path to gene expression data file
        pseudotime_file (str, optional): Path to pseudotime data file
        gene_names_file (str, optional): Path to gene names file
    
    Returns:
        tuple: (expression_data, pseudotime, gene_names)
    """
    # Load gene expression data
    try:
        if expression_file.endswith('.csv'):
            data = pd.read_csv(expression_file, index_col=0)
        elif expression_file.endswith('.tsv'):
            data = pd.read_csv(expression_file, sep='\t', index_col=0)
        else:
            # Try to load as CSV by default
            data = pd.read_csv(expression_file, index_col=0)
    except Exception as e:
        raise ValueError(f"Error loading expression data: {str(e)}")
    
    # Get gene names if available from data
    gene_names = data.index.tolist() if data.index.name == 'gene' else None
    
    # Convert to numpy array with genes as rows, cells as columns
    expression_data = data.values
    if gene_names is None:
        # If data is in cells x genes format, transpose it
        if data.shape[0] > data.shape[1]:
            expression_data = expression_data.T
    
    # Load pseudotime if file provided
    pseudotime = None
    if pseudotime_file:
        try:
            if pseudotime_file.endswith('.csv'):
                time_data = pd.read_csv(pseudotime_file)
            elif pseudotime_file.endswith('.tsv'):
                time_data = pd.read_csv(pseudotime_file, sep='\t')
            else:
                time_data = pd.read_csv(pseudotime_file)
            
            # Extract pseudotime column
            time_col = [col for col in time_data.columns if 'time' in col.lower() or 'pseudotime' in col.lower()]
            if time_col:
                pseudotime = time_data[time_col[0]].values
            else:
                # If no time column found, use the second column (assuming first is cell ID)
                pseudotime = time_data.iloc[:, 1].values
        except Exception as e:
            raise ValueError(f"Error loading pseudotime data: {str(e)}")
    
    # If pseudotime still not found, generate linear pseudotime
    if pseudotime is None or len(pseudotime) == 0:
        print("No pseudotime provided, generating linear pseudotime")
        pseudotime = np.linspace(0, 1, expression_data.shape[1])
    
    # Load gene names if file provided
    if gene_names is None and gene_names_file:
        try:
            with open(gene_names_file, 'r') as f:
                gene_names = [line.strip() for line in f.readlines()]
        except Exception as e:
            print(f"Warning: Could not load gene names: {str(e)}")
    
    # Create default gene names if not available
    if gene_names is None:
        gene_names = [f"Gene_{i}" for i in range(expression_data.shape[0])]
    
    print(f"Loaded data with {expression_data.shape[0]} genes and {expression_data.shape[1]} cells")
    return expression_data, pseudotime, gene_names

def run_pipeline(expression_data, pseudotime, gene_names=None, n_components=20, 
                time_param=1.0, n_measurements=1024, output_dir="output", 
                run_classical=True, top_n=20, run_trajectory=False):
    """
    Run the complete quantum differential gene expression analysis pipeline.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray): Pseudotime values for each cell
        gene_names (list, optional): List of gene names
        n_components (int): Number of principal components to use
        time_param (float): Time parameter for Hamiltonian evolution
        n_measurements (int): Number of quantum measurements to perform
        output_dir (str): Directory to save output files
        run_classical (bool): Whether to run classical analysis for comparison
        top_n (int): Number of top genes to compare between methods
        run_trajectory (bool): Whether to run trajectory analysis
    
    Returns:
        dict: Complete analysis results
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize results dictionary
    results = {}
    
    # Run quantum differential analysis
    print("\n==== Running Quantum Differential Analysis ====")
    start_time = time.time()
    quantum_results = quantum_differential_analysis(
        expression_data, pseudotime, n_components, time_param, n_measurements
    )
    end_time = time.time()
    quantum_results['execution_time'] = end_time - start_time
    print(f"Quantum analysis completed in {quantum_results['execution_time']:.2f} seconds")
    
    results['quantum'] = quantum_results
    
    # Run classical analysis for comparison if requested
    if run_classical:
        print("\n==== Running Classical Differential Analysis for Comparison ====")
        classical_results = classical_differential_analysis(
            expression_data, pseudotime, n_components
        )
        print(f"Classical analysis completed in {classical_results['execution_time']:.2f} seconds")
        
        # Compare methods
        print("\n==== Comparing Methods ====")
        comparison = compare_methods(quantum_results, classical_results, gene_names, top_n)
        print(f"Overlap between methods: {comparison['top_n_overlap_percentage']:.1f}% of top {top_n} genes")
        print(f"Rank correlation: {comparison['rank_correlation']:.3f}")
        
        results['classical'] = classical_results
        results['comparison'] = comparison
    
    # Run trajectory analysis if requested
    if run_trajectory:
        print("\n==== Running Quantum Trajectory Analysis ====")
        trajectory_output_dir = os.path.join(output_dir, "trajectory")
        os.makedirs(trajectory_output_dir, exist_ok=True)
        
        start_time = time.time()
        trajectory_results = quantum_trajectory_analysis(
            expression_data, 
            pseudotime, 
            gene_names,
            n_components=n_components,
            time_param=time_param,
            n_measurements=n_measurements,
            quantum_backend='qiskit',  # Default to qiskit
            output_dir=trajectory_output_dir
        )
        end_time = time.time()
        
        trajectory_execution_time = end_time - start_time
        print(f"Trajectory analysis completed in {trajectory_execution_time:.2f} seconds")
        
        results['trajectory'] = {
            'results': trajectory_results,
            'execution_time': trajectory_execution_time
        }
    
    # Save top differentially expressed genes
    save_top_genes(quantum_results, gene_names, output_dir, "quantum_top_genes.csv")
    if run_classical:
        save_top_genes({'diff_results': {'diff_genes_indices': classical_results['diff_genes_indices'], 
                                        'diff_scores': classical_results['combined_scores']}}, 
                     gene_names, output_dir, "classical_top_genes.csv")
    
    # Generate and save visualizations
    generate_visualizations(results, expression_data, pseudotime, gene_names, output_dir)
    
    return results

def save_top_genes(results, gene_names, output_dir, filename, top_n=50):
    """
    Save top differentially expressed genes to a CSV file.
    
    Args:
        results (dict): Analysis results
        gene_names (list): List of gene names
        output_dir (str): Directory to save output file
        filename (str): Output filename
        top_n (int): Number of top genes to save
    """
    # Get top genes and scores
    top_indices = results['diff_results']['diff_genes_indices'][:top_n]
    top_scores = results['diff_results']['diff_scores'][top_indices]
    
    # Create dataframe
    top_genes_df = pd.DataFrame({
        'Rank': range(1, len(top_indices) + 1),
        'Gene_Index': top_indices,
        'Gene_Name': [gene_names[i] for i in top_indices] if gene_names else top_indices,
        'Score': top_scores
    })
    
    # Save to CSV
    output_path = os.path.join(output_dir, filename)
    top_genes_df.to_csv(output_path, index=False)
    print(f"Saved top genes to {output_path}")

def generate_visualizations(results, expression_data, pseudotime, gene_names, output_dir):
    """
    Generate and save visualizations of analysis results.
    
    Args:
        results (dict): Complete analysis results
        expression_data (np.ndarray): Gene expression data
        pseudotime (np.ndarray): Pseudotime values
        gene_names (list): List of gene names
        output_dir (str): Directory to save visualizations
    """
    quantum_results = results['quantum']
    
    # 1. Plot top eigenvalues
    plt.figure(figsize=(10, 6))
    eigenvalues = quantum_results['quantum_signatures']['eigenvalues'][:10]  # Plot top 10
    plt.bar(range(len(eigenvalues)), eigenvalues)
    plt.xlabel('Index')
    plt.ylabel('Eigenvalue')
    plt.title('Top Eigenvalues from Quantum Analysis')
    plt.savefig(os.path.join(output_dir, 'eigenvalues.png'))
    plt.close()
    
    # 2. Plot probability distribution of quantum states
    top_states = quantum_results['quantum_signatures']['top_states'][:10]  # Plot top 10
    states = [s[0] for s in top_states]
    probs = [s[1] for s in top_states]
    plt.figure(figsize=(12, 6))
    plt.bar(states, probs)
    plt.xticks(rotation=45)
    plt.xlabel('Quantum State')
    plt.ylabel('Probability')
    plt.title('Top Quantum States by Probability')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'quantum_states.png'))
    plt.close()
    
    # 3. Plot top differentially expressed genes along pseudotime
    plt.figure(figsize=(12, 8))
    # Sort by pseudotime
    sort_idx = np.argsort(pseudotime)
    sorted_time = pseudotime[sort_idx]
    sorted_expr = expression_data[:, sort_idx]
    
    # Plot top 5 genes
    top_gene_indices = quantum_results['diff_results']['diff_genes_indices'][:5]
    for i, gene_idx in enumerate(top_gene_indices):
        gene_name = gene_names[gene_idx] if gene_names else f"Gene_{gene_idx}"
        plt.plot(sorted_time, sorted_expr[gene_idx], label=gene_name)
    
    plt.xlabel('Pseudotime')
    plt.ylabel('Expression Level')
    plt.title('Top Differentially Expressed Genes Along Pseudotime')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'top_genes_expression.png'))
    plt.close()
    
    # 4. If classical results available, compare execution times
    if 'classical' in results:
        plt.figure(figsize=(8, 6))
        methods = ['Quantum', 'Classical']
        times = [quantum_results['execution_time'], results['classical']['execution_time']]
        plt.bar(methods, times)
        plt.ylabel('Execution Time (seconds)')
        plt.title('Method Performance Comparison')
        for i, v in enumerate(times):
            plt.text(i, v + 0.1, f"{v:.2f}s", ha='center')
        plt.savefig(os.path.join(output_dir, 'execution_time_comparison.png'))
        plt.close()
    
    # 5. Plot gene weights from quantum analysis
    plt.figure(figsize=(10, 6))
    eigen_weights = quantum_results['diff_results']['eigen_weights']
    # Show distribution of weights
    plt.hist(eigen_weights, bins=30)
    plt.xlabel('Quantum Weight')
    plt.ylabel('Frequency')
    plt.title('Distribution of Gene Weights from Quantum Analysis')
    plt.savefig(os.path.join(output_dir, 'gene_weights_distribution.png'))
    plt.close()
    
    # Save figures as a collection for the report
    print(f"Saved visualizations to {output_dir}/")
