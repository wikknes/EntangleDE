import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from qiskit_aer import Aer
from qiskit.primitives import Sampler

def create_hamiltonian(gene_expression_data, pseudotime):
    """
    Create a Hamiltonian embedding for gene expression data along pseudotime.
    
    Args:
        gene_expression_data (np.ndarray): Matrix of gene expression values
                                          (genes x cells)
        pseudotime (np.ndarray): Vector of pseudotime values for each cell
    
    Returns:
        np.ndarray: Hamiltonian matrix representing the gene expression dynamics
    """
    n_genes, n_cells = gene_expression_data.shape
    
    # Sort data by pseudotime
    sorted_indices = np.argsort(pseudotime)
    sorted_data = gene_expression_data[:, sorted_indices]
    sorted_time = pseudotime[sorted_indices]
    
    # Calculate time differences
    dt = np.diff(sorted_time)
    dt = np.append(dt, dt[-1])  # Repeat last time difference for the last point
    
    # Calculate transition rates between consecutive time points
    transition_rates = np.zeros((n_genes, n_genes))
    
    for t in range(n_cells-1):
        # Calculate expression differences
        dE = (sorted_data[:, t+1] - sorted_data[:, t]) / dt[t]
        
        # Update transition rates based on expression changes
        for i in range(n_genes):
            for j in range(n_genes):
                # Simplified model: transition rate proportional to expression level
                # and expression change
                transition_rates[i, j] += sorted_data[j, t] * dE[i]
    
    # Normalize transition rates
    transition_rates = transition_rates / n_cells
    
    # Ensure Hermitian property for Hamiltonian
    hamiltonian = (transition_rates + transition_rates.T) / 2
    
    return hamiltonian

def encode_hamiltonian_to_circuit(hamiltonian, time_evolution=1.0, num_qubits=None):
    """
    Encode Hamiltonian into a quantum circuit for time evolution.
    
    Args:
        hamiltonian (np.ndarray): Hamiltonian matrix
        time_evolution (float): Time parameter for evolution
        num_qubits (int, optional): Number of qubits to use, defaults to log2(matrix_size)
    
    Returns:
        QuantumCircuit: Quantum circuit implementing Hamiltonian evolution
    """
    n_genes = hamiltonian.shape[0]
    
    # Determine number of qubits needed
    if num_qubits is None:
        num_qubits = int(np.ceil(np.log2(n_genes)))
    
    # Pad Hamiltonian if necessary
    padded_size = 2**num_qubits
    if padded_size > n_genes:
        padded_hamiltonian = np.zeros((padded_size, padded_size))
        padded_hamiltonian[:n_genes, :n_genes] = hamiltonian
        hamiltonian = padded_hamiltonian
    
    # Create quantum circuit
    circuit = QuantumCircuit(num_qubits)
    
    # Calculate time evolution operator - ensure it's unitary
    # Use eigendecomposition to ensure unitarity
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    # Create unitary matrix using eigendecomposition
    unitary_matrix = eigenvectors @ np.diag(np.exp(-1j * time_evolution * eigenvalues)) @ eigenvectors.conj().T
    
    # Convert to Operator
    evolution_operator = Operator(unitary_matrix)
    
    # Apply evolution operator to circuit
    circuit.unitary(evolution_operator, range(num_qubits), label='H_evolution')
    
    return circuit

def hamiltonian_embedding(gene_expression_data, pseudotime, time_param=1.0, num_qubits=None):
    """
    Perform Hamiltonian embedding of gene expression data along pseudotime.
    
    Args:
        gene_expression_data (np.ndarray): Matrix of gene expression values
                                          (genes x cells)
        pseudotime (np.ndarray): Vector of pseudotime values for each cell
        time_param (float): Time parameter for Hamiltonian evolution
        num_qubits (int, optional): Number of qubits to use, defaults to log2(n_genes)
    
    Returns:
        tuple: (hamiltonian, quantum_circuit)
    """
    # Create Hamiltonian from gene expression data
    hamiltonian = create_hamiltonian(gene_expression_data, pseudotime)
    
    # Encode Hamiltonian into quantum circuit
    circuit = encode_hamiltonian_to_circuit(hamiltonian, time_param, num_qubits)
    
    return hamiltonian, circuit
