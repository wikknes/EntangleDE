# EntangleDE: Quantum-Enhanced Gene Expression & Trajectory Analysis

## Overview

EntangleDE is a computational pipeline for analyzing gene expression data using quantum computing approaches. This innovative toolkit applies Hamiltonian embedding to model gene expression dynamics, providing enhanced capabilities for differential expression analysis and trajectory inference in single-cell RNA sequencing data.

## Key Features

- **Quantum Gene Expression Analysis**:
  - Quantum Hamiltonian embedding of gene expression data
  - Analysis of gene expression changes along pseudotime trajectories
  - Identification of differentially expressed genes with quantum advantage
  - Benchmarking against classical differential expression methods

- **Quantum Trajectory Analysis**:
  - Trajectory inference for time series scRNA-seq data
  - Quantum clustering for identifying cell states
  - Force-directed graph construction using quantum optimization
  - Cell ordering along developmental pathways
  - Pseudotime refinement with quantum algorithms
  
- **Visualization and Reporting**:
  - Comprehensive visualizations of analysis results
  - Gene expression dynamics along trajectories
  - Cell state transition mapping
  - Quality metrics for trajectory evaluation
  
- **User-Friendly Interface**:
  - Command-line interface for batch processing
  - Programmatic API for custom workflows
  - Integration with Scanpy/AnnData ecosystem

## Installation

```bash
# Clone the repository
git clone https://github.com/wikknes/EntangleDE.git
cd qdta

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### Differential Expression Analysis

```bash
# Run differential expression analysis
python main.py \
  --expression path/to/expression_data.csv \
  --pseudotime path/to/pseudotime.csv \
  --output results
```

### Trajectory Analysis

```bash
# Run trajectory analysis
python main.py \
  --expression path/to/expression_data.csv \
  --pseudotime path/to/pseudotime.csv \
  --output results \
  --trajectory
```

### Programmatic API

```python
from src.pipeline import load_data, run_pipeline

# Load your data
expression_data, pseudotime, gene_names = load_data(
    "path/to/expression_data.csv",
    "path/to/pseudotime.csv",
    "path/to/gene_names.txt"
)

# Run the analysis pipeline with trajectory inference
results = run_pipeline(
    expression_data, 
    pseudotime,
    gene_names,
    run_trajectory=True
)
```

## Documentation

- For detailed usage instructions and examples, see the [COOKBOOK.md](COOKBOOK.md) file.
- For trajectory analysis documentation, see the [TRAJECTORY_COOKBOOK.md](TRAJECTORY_COOKBOOK.md) file.
- For benchmarking information, see the [BENCHMARK.md](BENCHMARK.md) file.

## Requirements

- Python 3.7+
- Qiskit and qiskit-aer
- NumPy, Pandas, SciPy
- Scikit-learn, StatsModels
- Matplotlib, Seaborn
- Scanpy (for trajectory analysis)
- NetworkX (for graph algorithms)

Optional:
- D-Wave Ocean SDK (for quantum annealing)

## License

This project is licensed under the MIT License - see the LICENSE file for details.


## Contact

For questions or collaboration opportunities, please contact: [vignesh.kumar@csir.res.in]
