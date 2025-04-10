import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import Aer
from qiskit.visualization import plot_histogram
from sklearn.decomposition import PCA
from src.hamiltonian_embedding import hamiltonian_embedding

def normalize_expression_data(expression_data):
    """
    Normalize gene expression data.
    
    Args:
        expression_data (np.ndarray): Raw gene expression data (genes x cells)
    
    Returns:
        np.ndarray: Normalized gene expression data
    """
    # Log transform (adding small constant to avoid log(0))
    log_data = np.log1p(expression_data)
    
    # Min-max normalization to [0,1] range per gene
    norm_data = np.zeros_like(log_data)
    for i in range(log_data.shape[0]):
        gene_min = np.min(log_data[i, :])
        gene_max = np.max(log_data[i, :])
        if gene_max > gene_min:
            norm_data[i, :] = (log_data[i, :] - gene_min) / (gene_max - gene_min)
        else:
            norm_data[i, :] = 0.5  # Set to middle value if no variation
    
    return norm_data

def perform_dimensionality_reduction(expression_data, n_components=20):
    """
    Reduce dimensionality of gene expression data.
    
    Args:
        expression_data (np.ndarray): Normalized gene expression data (genes x cells)
        n_components (int): Number of components to keep
    
    Returns:
        np.ndarray: Reduced dimension gene expression data (n_components x cells)
    """
    # Transpose data for PCA (cells x genes)
    data_t = expression_data.T
    
    # Perform PCA
    pca = PCA(n_components=min(n_components, min(data_t.shape)))
    reduced_data_t = pca.fit_transform(data_t)
    
    # Transpose back to genes x cells
    reduced_data = reduced_data_t.T
    
    print(f"Explained variance ratio: {np.sum(pca.explained_variance_ratio_):.3f}")
    
    return reduced_data, pca

def calculate_quantum_signatures(hamiltonian, circuit, n_measurements=1024):
    """
    Calculate quantum signatures from Hamiltonian circuit.
    
    Args:
        hamiltonian (np.ndarray): Hamiltonian matrix
        circuit (QuantumCircuit): Quantum circuit implementing the Hamiltonian
        n_measurements (int): Number of measurements to perform
    
    Returns:
        dict: Quantum signatures extracted from the circuit execution
    """
    # Add measurement operations
    meas_circuit = circuit.copy()
    meas_circuit.measure_all()
    
    # Execute circuit
    simulator = Aer.get_backend('qasm_simulator')
    job = simulator.run(meas_circuit, shots=n_measurements)
    result = job.result()
    counts = result.get_counts(meas_circuit)
    
    # Calculate quantum signatures
    num_qubits = circuit.num_qubits
    
    # Eigenvalue estimation (simplified)
    eigenvalues = np.linalg.eigvals(hamiltonian)
    sorted_eigenvalues = np.sort(np.real(eigenvalues))[::-1]  # Sort in descending order
    
    # State probabilities
    probabilities = {}
    for state, count in counts.items():
        probabilities[state] = count / n_measurements
    
    # Calculate entropy of measurement outcomes
    entropy = 0
    for prob in probabilities.values():
        if prob > 0:
            entropy -= prob * np.log2(prob)
    
    # Extract top states by probability
    sorted_states = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
    top_states = sorted_states[:min(10, len(sorted_states))]
    
    return {
        'eigenvalues': sorted_eigenvalues,
        'top_eigenvalue': sorted_eigenvalues[0] if len(sorted_eigenvalues) > 0 else 0,
        'entropy': entropy,
        'top_states': top_states,
        'probabilities': probabilities
    }

def find_differentially_expressed_genes(original_data, pseudotime, quantum_signatures, 
                                      time_points=5, pca_components=None, pca=None):
    """
    Identify differentially expressed genes using quantum signatures.
    
    Args:
        original_data (np.ndarray): Original gene expression data (genes x cells)
        pseudotime (np.ndarray): Pseudotime values for each cell
        quantum_signatures (dict): Quantum signatures from circuit execution
        time_points (int): Number of time points to analyze
        pca_components (np.ndarray, optional): Reduced dimension data if PCA was applied
        pca (PCA, optional): Fitted PCA object if dimensionality reduction was applied
    
    Returns:
        dict: Information about differentially expressed genes
    """
    n_genes, n_cells = original_data.shape
    
    # Divide pseudotime into segments
    time_bins = np.linspace(pseudotime.min(), pseudotime.max(), time_points+1)
    
    # Calculate average expression in each time bin
    time_bin_expressions = []
    for i in range(time_points):
        bin_mask = (pseudotime >= time_bins[i]) & (pseudotime < time_bins[i+1])
        if np.sum(bin_mask) > 0:
            bin_avg = np.mean(original_data[:, bin_mask], axis=1)
        else:
            bin_avg = np.zeros(n_genes)
        time_bin_expressions.append(bin_avg)
    
    # Calculate expression differences between consecutive time points
    diff_expressions = []
    for i in range(1, len(time_bin_expressions)):
        diff = time_bin_expressions[i] - time_bin_expressions[i-1]
        diff_expressions.append(diff)
    
    # Use eigenvalues to weight genes based on contribution to dynamics
    eigen_weights = None
    if pca is not None and pca_components is not None:
        # Map eigenvalues back to gene space through PCA loading matrix
        eigenvalues = quantum_signatures['eigenvalues']
        if len(eigenvalues) >= pca.components_.shape[0]:
            # Use top eigenvalues corresponding to number of PCA components
            eigen_subset = eigenvalues[:pca.components_.shape[0]]
            # Weight PCA components by eigenvalues and transform back to gene space
            weighted_components = eigen_subset[:, np.newaxis] * pca.components_
            eigen_weights = np.sum(np.abs(weighted_components), axis=0)
            
            # Normalize weights
            if np.max(eigen_weights) > 0:
                eigen_weights = eigen_weights / np.max(eigen_weights)
    
    # If PCA wasn't used or mapping failed, use a simpler approach
    if eigen_weights is None:
        # Calculate variance of expression across time points as a simple weight
        time_point_matrix = np.vstack(time_bin_expressions)
        eigen_weights = np.var(time_point_matrix, axis=0)
        if np.max(eigen_weights) > 0:
            eigen_weights = eigen_weights / np.max(eigen_weights)
    
    # Combine time-based differences with quantum eigenvalue weights
    # to score genes for differential expression
    diff_scores = np.zeros(n_genes)
    for diff in diff_expressions:
        diff_scores += np.abs(diff)
    
    # Final score combines difference magnitude and eigenvalue weighting
    final_scores = diff_scores * eigen_weights
    
    # Sort genes by final score
    sorted_indices = np.argsort(final_scores)[::-1]
    
    return {
        'diff_genes_indices': sorted_indices,
        'diff_scores': final_scores,
        'time_bin_expressions': time_bin_expressions,
        'eigen_weights': eigen_weights
    }

def quantum_differential_analysis(expression_data, pseudotime, 
                                 n_components=20, time_param=1.0, n_measurements=1024):
    """
    Perform quantum differential gene expression analysis along pseudotime.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray): Pseudotime values for each cell
        n_components (int): Number of principal components to use
        time_param (float): Time parameter for Hamiltonian evolution
        n_measurements (int): Number of quantum measurements to perform
    
    Returns:
        dict: Results of differential gene expression analysis
    """
    # Save original data dimensions
    n_genes, n_cells = expression_data.shape
    print(f"Analyzing {n_genes} genes across {n_cells} cells")
    
    # Normalize data
    normalized_data = normalize_expression_data(expression_data)
    
    # Apply dimensionality reduction if n_genes is large
    pca = None
    reduced_data = None
    if n_genes > n_components:
        print(f"Reducing dimensions from {n_genes} to {n_components} using PCA")
        reduced_data, pca = perform_dimensionality_reduction(normalized_data, n_components)
        data_for_hamiltonian = reduced_data
    else:
        data_for_hamiltonian = normalized_data
    
    # Apply Hamiltonian embedding
    print("Applying Hamiltonian embedding...")
    hamiltonian, circuit = hamiltonian_embedding(data_for_hamiltonian, pseudotime, time_param)
    
    # Calculate quantum signatures
    print("Executing quantum circuit and extracting signatures...")
    quantum_signatures = calculate_quantum_signatures(hamiltonian, circuit, n_measurements)
    
    # Find differentially expressed genes
    print("Identifying differentially expressed genes...")
    diff_results = find_differentially_expressed_genes(
        expression_data, pseudotime, quantum_signatures, 5, reduced_data, pca
    )
    
    # Compile results
    results = {
        'hamiltonian': hamiltonian,
        'circuit': circuit,
        'quantum_signatures': quantum_signatures,
        'diff_results': diff_results,
        'n_genes': n_genes,
        'n_cells': n_cells,
        'n_components_used': data_for_hamiltonian.shape[0]
    }
    
    return results
