import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
import time
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import mannwhitneyu

def classical_differential_analysis(expression_data, pseudotime, n_components=20):
    """
    Perform classical differential gene expression analysis along pseudotime.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray): Pseudotime values for each cell
        n_components (int): Number of principal components to use
    
    Returns:
        dict: Results of differential gene expression analysis
    """
    n_genes, n_cells = expression_data.shape
    print(f"Analyzing {n_genes} genes across {n_cells} cells using classical methods")
    
    # Start timing
    start_time = time.time()
    
    # Log transform data (add small constant to avoid log(0))
    log_data = np.log1p(expression_data)
    
    # Normalize data
    norm_data = np.zeros_like(log_data)
    for i in range(log_data.shape[0]):
        gene_min = np.min(log_data[i, :])
        gene_max = np.max(log_data[i, :])
        if gene_max > gene_min:
            norm_data[i, :] = (log_data[i, :] - gene_min) / (gene_max - gene_min)
        else:
            norm_data[i, :] = 0.5
    
    # Apply PCA if needed
    pca_results = None
    if n_genes > n_components:
        print(f"Reducing dimensions from {n_genes} to {n_components} using PCA")
        pca = PCA(n_components=min(n_components, min(norm_data.T.shape)))
        pca_results = pca.fit_transform(norm_data.T).T
        print(f"Explained variance ratio: {np.sum(pca.explained_variance_ratio_):.3f}")
    
    # Method 1: Spearman correlation with pseudotime
    correlation_scores = np.zeros(n_genes)
    for i in range(n_genes):
        corr, _ = spearmanr(expression_data[i, :], pseudotime)
        correlation_scores[i] = abs(corr) if not np.isnan(corr) else 0
    
    # Method 2: Smoothed expression changes
    # Divide pseudotime into segments
    time_points = 5
    time_bins = np.linspace(pseudotime.min(), pseudotime.max(), time_points+1)
    
    # Calculate average expression in each time bin
    time_bin_expressions = []
    for i in range(time_points):
        bin_mask = (pseudotime >= time_bins[i]) & (pseudotime < time_bins[i+1])
        if np.sum(bin_mask) > 0:
            bin_avg = np.mean(expression_data[:, bin_mask], axis=1)
        else:
            bin_avg = np.zeros(n_genes)
        time_bin_expressions.append(bin_avg)
    
    # Calculate expression changes between consecutive time points
    change_scores = np.zeros(n_genes)
    for i in range(1, len(time_bin_expressions)):
        change = np.abs(time_bin_expressions[i] - time_bin_expressions[i-1])
        change_scores += change
    
    # Method 3: LOWESS smoothing and derivation
    lowess_scores = np.zeros(n_genes)
    for i in range(n_genes):
        # Sort by pseudotime
        sorted_indices = np.argsort(pseudotime)
        sorted_expr = expression_data[i, sorted_indices]
        sorted_time = pseudotime[sorted_indices]
        
        # Apply LOWESS smoothing
        try:
            smoothed = lowess(sorted_expr, sorted_time, frac=0.3, it=1, return_sorted=False)
            # Calculate absolute derivative
            derivative = np.abs(np.diff(smoothed))
            lowess_scores[i] = np.mean(derivative)
        except:
            lowess_scores[i] = 0
    
    # Method 4: Early vs Late comparison
    # Compare expression between early and late pseudotime
    early_late_scores = np.zeros(n_genes)
    median_time = np.median(pseudotime)
    early_mask = pseudotime < median_time
    late_mask = pseudotime >= median_time
    
    for i in range(n_genes):
        early_expr = expression_data[i, early_mask]
        late_expr = expression_data[i, late_mask]
        try:
            # Use Mann-Whitney U test to compare distributions
            u_stat, p_value = mannwhitneyu(early_expr, late_expr, alternative='two-sided')
            # Convert p-value to score (lower p-value = higher score)
            early_late_scores[i] = -np.log10(p_value) if p_value > 0 else 0
        except:
            early_late_scores[i] = 0
    
    # Normalize scores
    correlation_scores = correlation_scores / np.max(correlation_scores) if np.max(correlation_scores) > 0 else correlation_scores
    change_scores = change_scores / np.max(change_scores) if np.max(change_scores) > 0 else change_scores
    lowess_scores = lowess_scores / np.max(lowess_scores) if np.max(lowess_scores) > 0 else lowess_scores
    early_late_scores = early_late_scores / np.max(early_late_scores) if np.max(early_late_scores) > 0 else early_late_scores
    
    # Combine scores (equal weights)
    combined_scores = (correlation_scores + change_scores + lowess_scores + early_late_scores) / 4
    
    # Sort genes by combined score
    sorted_indices = np.argsort(combined_scores)[::-1]
    
    # End timing
    end_time = time.time()
    execution_time = end_time - start_time
    
    results = {
        'diff_genes_indices': sorted_indices,
        'correlation_scores': correlation_scores,
        'change_scores': change_scores,
        'lowess_scores': lowess_scores,
        'early_late_scores': early_late_scores,
        'combined_scores': combined_scores,
        'time_bin_expressions': time_bin_expressions,
        'execution_time': execution_time
    }
    
    return results

def compare_methods(quantum_results, classical_results, gene_names=None, top_n=20):
    """
    Compare quantum and classical differential expression analysis results.
    
    Args:
        quantum_results (dict): Results from quantum differential analysis
        classical_results (dict): Results from classical differential analysis
        gene_names (list, optional): List of gene names
        top_n (int): Number of top genes to compare
    
    Returns:
        dict: Comparison metrics between methods
    """
    # Get top differentially expressed genes from each method
    quantum_top_genes = quantum_results['diff_results']['diff_genes_indices'][:top_n]
    classical_top_genes = classical_results['diff_genes_indices'][:top_n]
    
    # Calculate overlap between methods
    overlap_genes = set(quantum_top_genes).intersection(set(classical_top_genes))
    overlap_percentage = len(overlap_genes) / top_n * 100
    
    # Calculate Spearman correlation between gene rankings
    # First create a mapping of gene index to rank for each method
    quantum_ranks = {gene_idx: rank for rank, gene_idx in enumerate(quantum_results['diff_results']['diff_genes_indices'])}
    classical_ranks = {gene_idx: rank for rank, gene_idx in enumerate(classical_results['diff_genes_indices'])}
    
    # Get ranks for all genes from both methods
    n_genes = len(quantum_results['diff_results']['diff_scores'])
    quantum_ranking = [quantum_ranks.get(i, n_genes) for i in range(n_genes)]
    classical_ranking = [classical_ranks.get(i, n_genes) for i in range(n_genes)]
    
    # Calculate rank correlation
    rank_correlation, p_value = spearmanr(quantum_ranking, classical_ranking)
    
    # Compile report
    report = {
        'quantum_execution_time': quantum_results.get('execution_time', 'Not measured'),
        'classical_execution_time': classical_results['execution_time'],
        'top_n_overlap_count': len(overlap_genes),
        'top_n_overlap_percentage': overlap_percentage,
        'rank_correlation': rank_correlation,
        'rank_correlation_p_value': p_value,
        'quantum_top_genes': quantum_top_genes.tolist(),
        'classical_top_genes': classical_top_genes.tolist(),
        'overlap_genes': list(overlap_genes)
    }
    
    # Add gene names if provided
    if gene_names is not None:
        report['quantum_top_gene_names'] = [gene_names[i] for i in quantum_top_genes]
        report['classical_top_gene_names'] = [gene_names[i] for i in classical_top_genes]
        report['overlap_gene_names'] = [gene_names[i] for i in overlap_genes]
    
    return report
