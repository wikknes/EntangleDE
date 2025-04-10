import unittest
import numpy as np
import pandas as pd
import os
import shutil
from tempfile import TemporaryDirectory

from src.hamiltonian_embedding import hamiltonian_embedding
from src.quantum_gene_analysis import quantum_differential_analysis
from src.classical_benchmark import classical_differential_analysis
from src.pipeline import load_data, run_pipeline

class TestQuantumDifferentialAnalysis(unittest.TestCase):
    
    def setUp(self):
        # Create simple test data
        self.n_genes = 10
        self.n_cells = 20
        np.random.seed(42)  # For reproducibility
        
        # Create expression data with some genes having pseudotime correlation
        self.pseudotime = np.linspace(0, 1, self.n_cells)
        self.expression_data = np.random.rand(self.n_genes, self.n_cells) * 0.5
        
        # Make first few genes correlated with pseudotime
        for i in range(3):
            self.expression_data[i, :] = i * self.pseudotime + np.random.rand(self.n_cells) * 0.2
        
        # Create temporary directory for test outputs
        self.temp_dir = TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        
        # Create gene names
        self.gene_names = [f"Gene_{i}" for i in range(self.n_genes)]
    
    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()
    
    def test_hamiltonian_embedding(self):
        # Test Hamiltonian embedding function
        hamiltonian, circuit = hamiltonian_embedding(self.expression_data, self.pseudotime)
        
        # Check output dimensions
        self.assertEqual(hamiltonian.shape, (self.n_genes, self.n_genes))
        self.assertGreaterEqual(circuit.num_qubits, int(np.ceil(np.log2(self.n_genes))))
    
    def test_quantum_analysis(self):
        # Test quantum differential analysis
        results = quantum_differential_analysis(self.expression_data, self.pseudotime, 5)
        
        # Check result structure
        self.assertIn('hamiltonian', results)
        self.assertIn('circuit', results)
        self.assertIn('quantum_signatures', results)
        self.assertIn('diff_results', results)
        
        # Check differential expression results
        diff_results = results['diff_results']
        self.assertEqual(len(diff_results['diff_genes_indices']), self.n_genes)
        self.assertEqual(len(diff_results['diff_scores']), self.n_genes)
        
        # Check if genes with pseudotime correlation are highly ranked
        top_genes = diff_results['diff_genes_indices'][:3]
        for i in range(3):
            self.assertIn(i, top_genes)
    
    def test_classical_analysis(self):
        # Test classical differential analysis
        results = classical_differential_analysis(self.expression_data, self.pseudotime)
        
        # Check result structure
        self.assertIn('diff_genes_indices', results)
        self.assertIn('combined_scores', results)
        self.assertIn('execution_time', results)
        
        # Check if genes with pseudotime correlation are highly ranked
        top_genes = results['diff_genes_indices'][:3]
        for i in range(3):
            self.assertIn(i, top_genes)
    
    def test_pipeline(self):
        # Create test data files
        expression_file = os.path.join(self.output_dir, 'test_expression.csv')
        pseudotime_file = os.path.join(self.output_dir, 'test_pseudotime.csv')
        genes_file = os.path.join(self.output_dir, 'test_genes.txt')
        
        # Save expression data to CSV
        pd.DataFrame(self.expression_data, index=self.gene_names).to_csv(expression_file)
        
        # Save pseudotime to CSV
        pd.DataFrame({'pseudotime': self.pseudotime}).to_csv(pseudotime_file, index=False)
        
        # Save gene names to text file
        with open(genes_file, 'w') as f:
            for gene in self.gene_names:
                f.write(f"{gene}\n")
        
        # Test data loading
        expr, time, genes = load_data(expression_file, pseudotime_file, genes_file)
        self.assertEqual(expr.shape, self.expression_data.shape)
        self.assertEqual(len(time), len(self.pseudotime))
        self.assertEqual(genes, self.gene_names)
        
        # Test pipeline execution with reduced parameters for testing speed
        results = run_pipeline(
            self.expression_data, self.pseudotime, self.gene_names,
            n_components=5, n_measurements=100, output_dir=self.output_dir
        )
        
        # Check output files
        self.assertTrue(os.path.exists(os.path.join(self.output_dir, 'quantum_top_genes.csv')))
        self.assertTrue(os.path.exists(os.path.join(self.output_dir, 'classical_top_genes.csv')))
        self.assertTrue(os.path.exists(os.path.join(self.output_dir, 'eigenvalues.png')))

if __name__ == '__main__':
    unittest.main()
