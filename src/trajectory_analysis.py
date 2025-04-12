import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from qiskit import QuantumCircuit
from qiskit_aer import Aer
try:
    # For newer Qiskit versions
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit_algorithms import QAOA
    QAOA_IMPORT = 'qiskit_algorithms'
except ImportError:
    # For older Qiskit versions
    from qiskit.algorithms.optimizers import COBYLA
    from qiskit.algorithms import QAOA
    QAOA_IMPORT = 'qiskit.algorithms'
from qiskit.circuit.library import ZZFeatureMap
from qiskit.quantum_info import Operator
import matplotlib.pyplot as plt
try:
    from dwave.system import DWaveSampler, EmbeddingComposite
    DWAVE_AVAILABLE = True
except ImportError:
    DWAVE_AVAILABLE = False
import dimod
import networkx as nx
from src.hamiltonian_embedding import create_hamiltonian, encode_hamiltonian_to_circuit

class QuantumTrajectoryAnalysis:
    """
    Quantum-enhanced single-cell RNA-seq trajectory analysis.
    
    This class implements trajectory inference for scRNA-seq data using quantum computing
    approaches, specifically focused on cells with time series information.
    """
    
    def __init__(self, n_components=20, time_param=1.0, n_measurements=1024, 
                 quantum_backend='qiskit', n_neighbors=10):
        """
        Initialize the trajectory analysis tool.
        
        Args:
            n_components (int): Number of principal components to use
            time_param (float): Time parameter for Hamiltonian evolution
            n_measurements (int): Number of quantum measurements to perform
            quantum_backend (str): 'qiskit' or 'dwave'
            n_neighbors (int): Number of nearest neighbors for graph construction
        """
        self.n_components = n_components
        self.time_param = time_param
        self.n_measurements = n_measurements
        self.quantum_backend = quantum_backend
        self.n_neighbors = n_neighbors
        
        # Results storage
        self.results = {}
        
    def preprocess_data(self, adata):
        """
        Preprocess AnnData object for trajectory analysis.
        
        Args:
            adata (AnnData): AnnData object containing gene expression data
            
        Returns:
            AnnData: Preprocessed AnnData object
        """
        print("Preprocessing data...")
        
        # Basic Scanpy preprocessing workflow
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(1000, adata.shape[1]))
        adata = adata[:, adata.var.highly_variable]
        
        # PCA for dimensionality reduction
        sc.pp.pca(adata, n_comps=self.n_components)
        
        # KNN computation for graph construction
        sc.pp.neighbors(adata, n_neighbors=self.n_neighbors)
        
        return adata
    
    def _construct_similarity_graph(self, adata, time_series=None):
        """
        Construct a similarity graph between cells based on gene expression and time information.
        
        Args:
            adata (AnnData): Preprocessed AnnData object
            time_series (np.ndarray, optional): Time points for each cell
            
        Returns:
            tuple: (similarity_matrix, graph)
        """
        # Get PCA coordinates
        X = adata.obsm['X_pca']
        
        # Calculate similarity based on PCA distance
        nbrs = NearestNeighbors(n_neighbors=self.n_neighbors).fit(X)
        distances, indices = nbrs.kneighbors(X)
        
        # Create similarity matrix (sparse)
        n_cells = X.shape[0]
        similarity = sparse.lil_matrix((n_cells, n_cells), dtype=float)
        
        # Fill similarity matrix with nearest neighbor connections
        for i in range(n_cells):
            similarity[i, indices[i]] = 1.0 - distances[i] / np.max(distances[i])
        
        # If time series is provided, adjust similarity by time proximity
        if time_series is not None:
            # Normalize time to [0,1]
            norm_time = (time_series - np.min(time_series)) / (np.max(time_series) - np.min(time_series))
            
            # Penalize connections between temporally distant cells
            for i in range(n_cells):
                for j in indices[i]:
                    time_diff = abs(norm_time[i] - norm_time[j])
                    time_penalty = np.exp(-5 * time_diff)  # Exponential penalty
                    similarity[i, j] *= time_penalty
        
        # Make symmetric
        similarity = (similarity + similarity.T) / 2
        
        # Convert to dense for quantum algorithms
        similarity_dense = similarity.toarray()
        
        # Create networkx graph
        graph = nx.from_numpy_array(similarity_dense)
        
        return similarity_dense, graph
    
    def _formulate_qubo(self, similarity_matrix, n_clusters=5):
        """
        Formulate the clustering problem as a QUBO for quantum optimization.
        
        Args:
            similarity_matrix (np.ndarray): Cell similarity matrix
            n_clusters (int): Number of clusters
            
        Returns:
            dict: QUBO dictionary
        """
        n_cells = similarity_matrix.shape[0]
        
        # Initialize QUBO dictionary
        Q = {}
        
        # Add terms that maximize similarity within clusters
        for i in range(n_cells):
            for j in range(i, n_cells):
                for k in range(n_clusters):
                    # Variable indices for cell i and j being in cluster k
                    idx_i_k = (i, k)
                    idx_j_k = (j, k)
                    
                    # Encourage cells with high similarity to be in the same cluster
                    # For each pair of cells, add a negative term to promote grouping
                    # similar cells together in the same cluster
                    if i == j:
                        # Diagonal terms
                        Q[(idx_i_k, idx_i_k)] = Q.get((idx_i_k, idx_i_k), 0.0) - 1.0
                    else:
                        # Off-diagonal terms: similarity encourages same cluster assignment
                        Q[(idx_i_k, idx_j_k)] = Q.get((idx_i_k, idx_j_k), 0.0) - similarity_matrix[i, j]
        
        # Add constraint: each cell belongs to exactly one cluster
        lambda_constraint = 5.0  # Penalty strength
        for i in range(n_cells):
            # Add quadratic terms to penalize multiple cluster assignments
            for k1 in range(n_clusters):
                for k2 in range(k1+1, n_clusters):
                    idx_i_k1 = (i, k1)
                    idx_i_k2 = (i, k2)
                    Q[(idx_i_k1, idx_i_k2)] = Q.get((idx_i_k1, idx_i_k2), 0.0) + lambda_constraint
            
            # Add linear terms to encourage exactly one cluster assignment
            for k in range(n_clusters):
                idx_i_k = (i, k)
                Q[(idx_i_k, idx_i_k)] = Q.get((idx_i_k, idx_i_k), 0.0) + lambda_constraint * (1 - 2)
        
        return Q
    
    def _solve_with_quantum_annealing(self, qubo, n_cells, n_clusters):
        """
        Solve the QUBO problem using D-Wave quantum annealer.
        
        Args:
            qubo (dict): QUBO dictionary
            n_cells (int): Number of cells
            n_clusters (int): Number of clusters
            
        Returns:
            np.ndarray: Cluster assignments for each cell
        """
        print("Solving with D-Wave quantum annealer...")
        
        try:
            if DWAVE_AVAILABLE:
                # Create quantum sampler
                sampler = EmbeddingComposite(DWaveSampler())
                
                # Run quantum annealing
                response = sampler.sample_qubo(qubo, num_reads=100)
                solution = response.first.sample
            else:
                raise ImportError("D-Wave Ocean SDK not available")
            
            # Extract cluster assignments
            cluster_assignments = np.zeros(n_cells, dtype=int)
            for (cell, cluster), value in solution.items():
                if value == 1:
                    cluster_assignments[cell] = cluster
            
            return cluster_assignments
            
        except Exception as e:
            print(f"D-Wave solver error: {e}")
            print("Falling back to classical simulated annealing")
            
            # Fallback to simulated annealing
            sampler = dimod.SimulatedAnnealingSampler()
            response = sampler.sample_qubo(qubo, num_reads=100)
            solution = response.first.sample
            
            # Extract cluster assignments
            cluster_assignments = np.zeros(n_cells, dtype=int)
            for (cell, cluster), value in solution.items():
                if value == 1:
                    cluster_assignments[cell] = cluster
            
            return cluster_assignments
    
    def _solve_with_qaoa(self, similarity_matrix, n_clusters):
        """
        Solve the clustering problem using Quantum Approximate Optimization Algorithm.
        
        Args:
            similarity_matrix (np.ndarray): Cell similarity matrix
            n_clusters (int): Number of clusters
            
        Returns:
            np.ndarray: Cluster assignments for each cell
        """
        print("Solving with QAOA...")
        
        try:
            n_cells = similarity_matrix.shape[0]
            
            # Import appropriate modules based on Qiskit version
            try:
                # For newer Qiskit versions
                from qiskit_algorithms.optimizers import COBYLA
                from qiskit_algorithms import QAOA
                from qiskit.primitives import Sampler
                from qiskit_optimization import QuadraticProgram
                from qiskit_optimization.converters import QuadraticProgramToQubo
                
                # Create a quadratic program
                qp = QuadraticProgram()
                
                # Add binary variables for cell-cluster assignments
                for i in range(n_cells):
                    for k in range(n_clusters):
                        qp.binary_var(name=f'x_{i}_{k}')
                
                # Set the objective (maximize similarity within clusters)
                linear = {}
                quadratic = {}
                
                # Add terms that maximize similarity within clusters
                for i in range(n_cells):
                    for j in range(i, n_cells):
                        for k in range(n_clusters):
                            # Variable names
                            var_i_k = f'x_{i}_{k}'
                            var_j_k = f'x_{j}_{k}'
                            
                            if i == j:
                                # Diagonal terms (single-cell clustering preference)
                                linear[var_i_k] = linear.get(var_i_k, 0.0) - 1.0
                            else:
                                # Off-diagonal terms (pairwise cell similarity)
                                quadratic[(var_i_k, var_j_k)] = quadratic.get((var_i_k, var_j_k), 0.0) - similarity_matrix[i, j]
                
                # Add constraints: each cell belongs to exactly one cluster
                for i in range(n_cells):
                    qp.linear_constraint(
                        linear={f'x_{i}_{k}': 1 for k in range(n_clusters)},
                        sense='==',
                        rhs=1,
                        name=f'cell_{i}_constraint'
                    )
                
                # Convert to QUBO
                qp2qubo = QuadraticProgramToQubo()
                qubo = qp2qubo.convert(qp)
                
                # Set up QAOA
                sampler = Sampler()
                optimizer = COBYLA(maxiter=100)
                qaoa = QAOA(sampler=sampler, optimizer=optimizer, reps=2)
                
                # Solve
                result = qaoa.compute_minimum_eigenvalue(qubo.to_ising()[0])
                
            except (ImportError, AttributeError):
                # For older Qiskit versions
                from qiskit.algorithms import QAOA
                from qiskit.algorithms.optimizers import COBYLA
                from qiskit.utils import QuantumInstance
                from qiskit_optimization import QuadraticProgram
                from qiskit_optimization.converters import QuadraticProgramToQubo
                
                # Create a quadratic program
                qp = QuadraticProgram()
                
                # Add binary variables for cell-cluster assignments
                for i in range(n_cells):
                    for k in range(n_clusters):
                        qp.binary_var(name=f'x_{i}_{k}')
                
                # Set the objective (maximize similarity within clusters)
                linear = {}
                quadratic = {}
                
                # Add terms that maximize similarity within clusters
                for i in range(n_cells):
                    for j in range(i, n_cells):
                        for k in range(n_clusters):
                            # Variable names
                            var_i_k = f'x_{i}_{k}'
                            var_j_k = f'x_{j}_{k}'
                            
                            if i == j:
                                # Diagonal terms (single-cell clustering preference)
                                linear[var_i_k] = linear.get(var_i_k, 0.0) - 1.0
                            else:
                                # Off-diagonal terms (pairwise cell similarity)
                                quadratic[(var_i_k, var_j_k)] = quadratic.get((var_i_k, var_j_k), 0.0) - similarity_matrix[i, j]
                
                # Add constraints: each cell belongs to exactly one cluster
                for i in range(n_cells):
                    qp.linear_constraint(
                        linear={f'x_{i}_{k}': 1 for k in range(n_clusters)},
                        sense='==',
                        rhs=1,
                        name=f'cell_{i}_constraint'
                    )
                
                # Convert to QUBO
                qp2qubo = QuadraticProgramToQubo()
                qubo = qp2qubo.convert(qp)
                
                # Set up QAOA
                quantum_instance = QuantumInstance(Aer.get_backend('qasm_simulator'), shots=1024)
                optimizer = COBYLA(maxiter=100)
                qaoa = QAOA(optimizer=optimizer, quantum_instance=quantum_instance, reps=2)
                
                # Solve
                result = qaoa.compute_minimum_eigenvalue(qubo.to_ising()[0])
            
            # Extract solution
            x = result.x
            
            # Convert to cluster assignments
            cluster_assignments = np.zeros(n_cells, dtype=int)
            for i in range(n_cells):
                for k in range(n_clusters):
                    var_idx = i * n_clusters + k
                    if x[var_idx] == 1:
                        cluster_assignments[i] = k
            
            return cluster_assignments
            
        except Exception as e:
            print(f"QAOA solver error: {e}")
            print("Falling back to classical clustering")
            
            # Fallback to classical spectral clustering
            from sklearn.cluster import SpectralClustering
            clustering = SpectralClustering(n_clusters=n_clusters, 
                                           affinity='precomputed',
                                           random_state=42)
            cluster_assignments = clustering.fit_predict(similarity_matrix)
            
            return cluster_assignments
    
    def _hamiltonian_based_trajectory(self, adata, time_series=None):
        """
        Use Hamiltonian evolution to infer trajectories.
        
        Args:
            adata (AnnData): Preprocessed AnnData object
            time_series (np.ndarray, optional): Time points for each cell
            
        Returns:
            tuple: (force_graph, pseudotime)
        """
        print("Performing Hamiltonian-based trajectory inference...")
        
        # Get PCA matrix (cells x features)
        X = adata.obsm['X_pca']
        
        # Transpose to get features x cells
        X_T = X.T
        
        # If time series not provided, initialize with PCA first component
        if time_series is None:
            time_series = X[:, 0]
            time_series = (time_series - np.min(time_series)) / (np.max(time_series) - np.min(time_series))
        
        # Create Hamiltonian from expression data and pseudotime
        hamiltonian = create_hamiltonian(X_T, time_series)
        
        # Encode into quantum circuit
        circuit = encode_hamiltonian_to_circuit(hamiltonian, self.time_param)
        
        # Execute circuit
        meas_circuit = circuit.copy()
        meas_circuit.measure_all()
        simulator = Aer.get_backend('qasm_simulator')
        job = simulator.run(meas_circuit, shots=self.n_measurements)
        result = job.result()
        counts = result.get_counts(meas_circuit)
        
        # Use Hamiltonian to create force-directed graph
        force_graph = nx.Graph()
        
        # Add nodes (cells)
        for i in range(X.shape[0]):
            force_graph.add_node(i, embedding=X[i])
        
        # Add edges with weights based on Hamiltonian
        for i in range(X.shape[0]):
            for j in range(i+1, X.shape[0]):
                # Edge weight inversely proportional to Hamiltonian element
                if i < hamiltonian.shape[0] and j < hamiltonian.shape[0]:
                    # Use absolute value of Hamiltonian element
                    h_ij = abs(hamiltonian[i, j])
                    if h_ij > 1e-5:  # Only add significant connections
                        force_graph.add_edge(i, j, weight=1.0/h_ij)
        
        # Calculate refined pseudotime based on force graph
        # Use shortest path from earliest cell (by original time) to all others
        sorted_indices = np.argsort(time_series)
        start_node = sorted_indices[0]  # Earliest cell
        
        # Use shortest path lengths as pseudotime
        path_lengths = nx.single_source_shortest_path_length(force_graph, start_node)
        
        # Convert to numpy array in original cell order
        refined_pseudotime = np.zeros(X.shape[0])
        for node, length in path_lengths.items():
            refined_pseudotime[node] = length
        
        # Normalize to [0,1]
        if np.max(refined_pseudotime) > 0:
            refined_pseudotime = refined_pseudotime / np.max(refined_pseudotime)
        
        return force_graph, refined_pseudotime
    
    def run_trajectory_analysis(self, adata, time_series=None, n_clusters=5):
        """
        Run trajectory analysis on the given AnnData object.
        
        Args:
            adata (AnnData): AnnData object containing gene expression data
            time_series (np.ndarray, optional): Time points for each cell
            n_clusters (int): Number of clusters for intermediate clustering
            
        Returns:
            dict: Analysis results
        """
        # Preprocess data
        adata = self.preprocess_data(adata)
        
        # Construct similarity graph between cells
        similarity_matrix, graph = self._construct_similarity_graph(adata, time_series)
        
        # Run quantum clustering depending on backend
        if self.quantum_backend == 'dwave':
            # Formulate as QUBO
            qubo = self._formulate_qubo(similarity_matrix, n_clusters)
            # Solve with quantum annealing
            cluster_assignments = self._solve_with_quantum_annealing(qubo, adata.shape[0], n_clusters)
        else:  # qiskit
            # Solve with QAOA
            cluster_assignments = self._solve_with_qaoa(similarity_matrix, n_clusters)
        
        # Add cluster assignments to AnnData
        adata.obs['quantum_clusters'] = cluster_assignments
        
        # Run trajectory inference using Hamiltonian approach
        force_graph, refined_pseudotime = self._hamiltonian_based_trajectory(adata, time_series)
        
        # Add refined pseudotime to AnnData
        adata.obs['quantum_pseudotime'] = refined_pseudotime
        
        # Store results
        self.results = {
            'adata': adata,
            'similarity_matrix': similarity_matrix,
            'graph': graph,
            'force_graph': force_graph,
            'quantum_clusters': cluster_assignments,
            'quantum_pseudotime': refined_pseudotime
        }
        
        return self.results
    
    def compute_trajectory_metrics(self, adata, true_pseudotime=None):
        """
        Compute evaluation metrics for the inferred trajectory.
        
        Args:
            adata (AnnData): AnnData object with inferred trajectory
            true_pseudotime (np.ndarray, optional): Ground truth pseudotime if available
            
        Returns:
            dict: Evaluation metrics
        """
        metrics = {}
        
        # Connectivity score (based on graph connectivity)
        connectivity = nx.average_node_connectivity(self.results['force_graph'])
        metrics['connectivity'] = connectivity
        
        # Kendall's tau correlation with true pseudotime if available
        if true_pseudotime is not None:
            from scipy.stats import kendalltau
            tau, p_value = kendalltau(self.results['quantum_pseudotime'], true_pseudotime)
            metrics['kendall_tau'] = tau
            metrics['kendall_p_value'] = p_value
        
        # Pseudo-time order stability (bootstrap test on subsamples)
        stability_scores = []
        n_bootstrap = 10
        
        X = adata.obsm['X_pca']
        n_cells = X.shape[0]
        
        for _ in range(n_bootstrap):
            # Subsample 80% of cells
            subsample_indices = np.random.choice(n_cells, size=int(0.8 * n_cells), replace=False)
            
            # Construct subsampled similarity matrix
            X_sub = X[subsample_indices]
            nbrs = NearestNeighbors(n_neighbors=min(self.n_neighbors, len(subsample_indices)-1)).fit(X_sub)
            distances, indices = nbrs.kneighbors(X_sub)
            
            sub_similarity = np.zeros((len(subsample_indices), len(subsample_indices)))
            for i in range(len(subsample_indices)):
                sub_similarity[i, indices[i]] = 1.0 - distances[i] / np.max(distances[i])
            
            # Make symmetric
            sub_similarity = (sub_similarity + sub_similarity.T) / 2
            
            # Create subsampled graph
            sub_graph = nx.from_numpy_array(sub_similarity)
            
            # Calculate pseudotime on subsampled graph
            start_node = 0  # Arbitrarily use first node as start
            path_lengths = nx.single_source_shortest_path_length(sub_graph, start_node)
            
            # Convert to numpy array
            sub_pseudotime = np.zeros(len(subsample_indices))
            for node, length in path_lengths.items():
                sub_pseudotime[node] = length
            
            # Normalize
            if np.max(sub_pseudotime) > 0:
                sub_pseudotime = sub_pseudotime / np.max(sub_pseudotime)
            
            # Compare to original pseudotime for these cells
            original_pt = self.results['quantum_pseudotime'][subsample_indices]
            # Calculate correlation
            from scipy.stats import spearmanr
            rho, _ = spearmanr(sub_pseudotime, original_pt)
            stability_scores.append(rho)
        
        # Average stability score
        metrics['stability'] = np.mean(stability_scores)
        
        return metrics
    
    def plot_trajectory(self, adata, embedding_key='X_pca', save_path=None):
        """
        Visualize the inferred trajectory.
        
        Args:
            adata (AnnData): AnnData object with inferred trajectory
            embedding_key (str): Key for cell embedding to use for visualization
            save_path (str, optional): Path to save the figure
            
        Returns:
            plt.Figure: Matplotlib figure object
        """
        import scanpy as sc
        
        if 'quantum_pseudotime' not in adata.obs:
            raise ValueError("Run trajectory analysis first")
        
        # Create UMAP if not already present
        if 'X_umap' not in adata.obsm:
            sc.pp.neighbors(adata, use_rep=embedding_key)
            sc.tl.umap(adata)
        
        # Get embedding
        if embedding_key == 'X_umap' and 'X_umap' in adata.obsm:
            embedding = adata.obsm['X_umap']
        else:
            embedding = adata.obsm[embedding_key][:, :2]  # Use first two components
        
        # Create figure
        fig, axs = plt.subplots(1, 2, figsize=(16, 7))
        
        # Plot clusters
        axs[0].scatter(embedding[:, 0], embedding[:, 1], 
                     c=adata.obs['quantum_clusters'], 
                     cmap='tab10', s=50, alpha=0.7)
        axs[0].set_title('Quantum Clustering')
        axs[0].set_xlabel(f'{embedding_key.replace("_", " ")} 1')
        axs[0].set_ylabel(f'{embedding_key.replace("_", " ")} 2')
        
        # Plot pseudotime
        scatter = axs[1].scatter(embedding[:, 0], embedding[:, 1], 
                          c=adata.obs['quantum_pseudotime'], 
                          cmap='viridis', s=50, alpha=0.7)
        axs[1].set_title('Quantum Pseudotime')
        axs[1].set_xlabel(f'{embedding_key.replace("_", " ")} 1')
        axs[1].set_ylabel(f'{embedding_key.replace("_", " ")} 2')
        plt.colorbar(scatter, ax=axs[1], label='Pseudotime')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_gene_dynamics(self, adata, genes, save_path=None):
        """
        Plot gene expression dynamics along the inferred trajectory.
        
        Args:
            adata (AnnData): AnnData object with inferred trajectory
            genes (list): List of gene names/indices to plot
            save_path (str, optional): Path to save the figure
            
        Returns:
            plt.Figure: Matplotlib figure object
        """
        from scipy import sparse
        from statsmodels.nonparametric.smoothers_lowess import lowess
        
        if 'quantum_pseudotime' not in adata.obs:
            raise ValueError("Run trajectory analysis first")
        
        # Get pseudotime
        pseudotime = adata.obs['quantum_pseudotime'].values
        
        # Sort cells by pseudotime
        sorted_indices = np.argsort(pseudotime)
        sorted_time = pseudotime[sorted_indices]
        
        # Create figure
        n_genes = len(genes)
        n_cols = min(3, n_genes)
        n_rows = int(np.ceil(n_genes / n_cols))
        
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
        if n_genes == 1:
            axs = np.array([axs])
        axs = axs.flatten()
        
        # Plot each gene
        for i, gene in enumerate(genes):
            if i < len(axs):
                # Get gene expression
                if isinstance(gene, int):
                    gene_expr = adata.X[:, gene]
                    gene_name = f"Gene {gene}"
                else:
                    if gene in adata.var_names:
                        gene_idx = adata.var_names.get_loc(gene)
                        gene_expr = adata.X[:, gene_idx]
                        gene_name = gene
                    else:
                        axs[i].text(0.5, 0.5, f"Gene {gene} not found", 
                                   horizontalalignment='center',
                                   verticalalignment='center')
                        continue
                
                # Convert to dense if sparse
                if sparse.issparse(adata.X):
                    gene_expr = gene_expr.toarray().flatten()
                
                # Scatter plot
                axs[i].scatter(pseudotime, gene_expr, alpha=0.5)
                
                # LOWESS smoothing
                try:
                    z = lowess(gene_expr, pseudotime, frac=0.3, it=1)
                    axs[i].plot(z[:, 0], z[:, 1], 'r-', linewidth=2)
                except Exception as e:
                    print(f"Warning: Could not apply LOWESS smoothing for {gene_name}: {e}")
                
                axs[i].set_title(gene_name)
                axs[i].set_xlabel('Pseudotime')
                axs[i].set_ylabel('Expression')
        
        # Hide unused subplots
        for i in range(n_genes, len(axs)):
            axs[i].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def export_trajectory(self, output_path):
        """
        Export trajectory analysis results to file.
        
        Args:
            output_path (str): Path to save the results
            
        Returns:
            None
        """
        if not self.results:
            raise ValueError("Run trajectory analysis first")
        
        # Export AnnData
        adata = self.results['adata']
        adata.write(f"{output_path}/trajectory_adata.h5ad")
        
        # Export graph as node-edge list
        graph = self.results['force_graph']
        nx.write_weighted_edgelist(graph, f"{output_path}/trajectory_graph.edgelist")
        
        # Export clusters and pseudotime as CSV
        pd.DataFrame({
            'Cell_ID': adata.obs_names,
            'Quantum_Cluster': self.results['quantum_clusters'],
            'Quantum_Pseudotime': self.results['quantum_pseudotime']
        }).to_csv(f"{output_path}/trajectory_results.csv", index=False)
        
        print(f"Trajectory analysis results exported to {output_path}")


def load_as_anndata(expression_data, pseudotime=None, gene_names=None, cell_ids=None):
    """
    Convert expression data to AnnData format for trajectory analysis.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray, optional): Pseudotime values for each cell
        gene_names (list, optional): List of gene names
        cell_ids (list, optional): List of cell identifiers
    
    Returns:
        AnnData: AnnData object ready for trajectory analysis
    """
    import scanpy as sc
    
    # Transpose to cells x genes format for AnnData
    X = expression_data.T
    
    # Create AnnData object
    adata = sc.AnnData(X)
    
    # Add gene names if provided
    if gene_names is not None:
        adata.var_names = gene_names
    
    # Add cell IDs if provided
    if cell_ids is not None:
        adata.obs_names = cell_ids
    else:
        adata.obs_names = [f"Cell_{i}" for i in range(X.shape[0])]
    
    # Add pseudotime if provided
    if pseudotime is not None:
        adata.obs['pseudotime'] = pseudotime
    
    return adata


def quantum_trajectory_analysis(expression_data, pseudotime=None, gene_names=None, 
                               n_components=20, time_param=1.0, n_measurements=1024,
                               quantum_backend='qiskit', n_clusters=5,
                               output_dir="trajectory_output"):
    """
    Perform quantum-enhanced trajectory analysis on gene expression data.
    
    Args:
        expression_data (np.ndarray): Gene expression data (genes x cells)
        pseudotime (np.ndarray, optional): Pseudotime values for each cell
        gene_names (list, optional): List of gene names
        n_components (int): Number of principal components to use
        time_param (float): Time parameter for Hamiltonian evolution
        n_measurements (int): Number of quantum measurements to perform
        quantum_backend (str): 'qiskit' or 'dwave'
        n_clusters (int): Number of clusters for intermediate clustering
        output_dir (str): Directory to save results
    
    Returns:
        dict: Trajectory analysis results
    """
    import os
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert to AnnData format
    adata = load_as_anndata(expression_data, pseudotime, gene_names)
    
    # Initialize trajectory analyzer
    analyzer = QuantumTrajectoryAnalysis(
        n_components=n_components,
        time_param=time_param,
        n_measurements=n_measurements,
        quantum_backend=quantum_backend
    )
    
    # Run trajectory analysis
    results = analyzer.run_trajectory_analysis(adata, pseudotime, n_clusters)
    
    # Compute evaluation metrics
    metrics = analyzer.compute_trajectory_metrics(results['adata'], pseudotime)
    
    # Create visualizations
    analyzer.plot_trajectory(results['adata'], save_path=f"{output_dir}/trajectory_plot.png")
    
    # If gene names are available, plot top 6 most variable genes
    if gene_names is not None:
        # Find top variable genes
        from scanpy import preprocessing as pp
        pp.highly_variable_genes(results['adata'])
        if np.sum(results['adata'].var['highly_variable']) > 0:
            top_genes = list(results['adata'].var_names[results['adata'].var['highly_variable']])[:6]
            analyzer.plot_gene_dynamics(results['adata'], top_genes, 
                                      save_path=f"{output_dir}/gene_dynamics.png")
    
    # Export results
    analyzer.export_trajectory(output_dir)
    
    # Return combined results
    return {
        'trajectory_results': results,
        'evaluation_metrics': metrics
    }