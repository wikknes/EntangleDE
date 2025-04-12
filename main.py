import argparse
import numpy as np
import pandas as pd
import os
from src.pipeline import load_data, run_pipeline

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='EntangleDE: Quantum-Enhanced Gene Expression & Trajectory Analysis')
    parser.add_argument('-e', '--expression', required=True, help='Path to gene expression data file')
    parser.add_argument('-p', '--pseudotime', help='Path to pseudotime data file')
    parser.add_argument('-g', '--genes', help='Path to gene names file')
    parser.add_argument('-c', '--components', type=int, default=20, help='Number of components for dimensionality reduction')
    parser.add_argument('-t', '--time', type=float, default=1.0, help='Time parameter for Hamiltonian evolution')
    parser.add_argument('-m', '--measurements', type=int, default=1024, help='Number of quantum measurements')
    parser.add_argument('-o', '--output', default='output', help='Output directory')
    parser.add_argument('--no-classical', action='store_true', help='Skip classical benchmark analysis')
    parser.add_argument('--top-n', type=int, default=20, help='Number of top genes to compare')
    
    # New trajectory analysis options
    parser.add_argument('--trajectory', action='store_true', help='Run trajectory analysis')
    parser.add_argument('--quantum-backend', choices=['qiskit', 'dwave'], default='qiskit', 
                      help='Quantum backend for trajectory analysis (default: qiskit)')
    parser.add_argument('--clusters', type=int, default=5, 
                      help='Number of clusters for trajectory analysis (default: 5)')
    
    args = parser.parse_args()
    
    # Load data
    expression_data, pseudotime, gene_names = load_data(
        args.expression, args.pseudotime, args.genes
    )
    
    # Run analysis pipeline
    results = run_pipeline(
        expression_data, 
        pseudotime, 
        gene_names,
        n_components=args.components,
        time_param=args.time,
        n_measurements=args.measurements,
        output_dir=args.output,
        run_classical=not args.no_classical,
        top_n=args.top_n,
        run_trajectory=args.trajectory
    )
    
    print(f"\nAnalysis complete! Results saved to {args.output}/")
    print(f"Top differentially expressed genes saved to {args.output}/quantum_top_genes.csv")
    
    if args.trajectory:
        print(f"Trajectory analysis results saved to {args.output}/trajectory/")

if __name__ == "__main__":
    main()
