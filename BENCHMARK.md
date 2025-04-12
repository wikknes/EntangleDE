# EntangleDE Performance Benchmarks

This document presents comprehensive performance benchmarks comparing the Quantum Differential Gene Expression Analysis (EntangleDE) method with classical approaches for identifying differentially expressed genes along pseudotime trajectories.

## Overview

We evaluated both approaches based on:
1. Execution time
2. Accuracy in identifying known differentially expressed genes
3. Overlap in top differentially expressed genes
4. Memory usage
5. Gene expression pattern detection

## Test Datasets

### Synthetic Datasets

We generated synthetic datasets with known differential expression patterns using the `generate_synthetic_data` function, which creates three types of patterns:

| Dataset | Genes | Cells | Diff. Genes | Expression Patterns |
|---------|-------|-------|-------------|---------------------|
| Small   | 20    | 50    | 5           | Linear, sigmoidal, bell-curve |
| Medium  | 100   | 200   | 10          | Linear, sigmoidal, bell-curve |
| Large   | 500   | 500   | 30          | Linear, sigmoidal, bell-curve |

The synthetic patterns include:
- **Linear**: Genes with expression that increases linearly with pseudotime
- **Sigmoidal**: Genes with an S-shaped expression curve along pseudotime
- **Bell-curve**: Genes with peak expression in the middle of the pseudotime trajectory

### Replication Instructions

To replicate these benchmarks, use the benchmark script:

```bash
# Small dataset
python benchmark.py --size small

# Medium dataset
python benchmark.py --size medium

# Large dataset
python benchmark.py --size large
```

## Latest Benchmark Results (2023)

The following results were obtained on the latest test run:

## Execution Time Comparison

Execution times measured on a MacBook:

| Dataset Size | Quantum Runtime | Classical Runtime | Speedup Factor |
|--------------|----------------|-------------------|----------------|
| Small (20 genes) | 0.02s | 0.12s | 6x |
| Medium (100 genes) | 0.05s | 0.57s | 11.4x |
| Large (500 genes) | 0.16s | 6.78s | 42.4x |

The quantum approach demonstrates remarkable performance advantages that increase with dataset size. For the large dataset, the quantum method is over 42 times faster than the classical approach, showing excellent scalability for larger problems.

> *Note: Quantum runtimes are from simulation on classical hardware using Qiskit. The actual performance will vary depending on hardware specifications.*

## Accuracy Assessment

Both methods were evaluated on their ability to identify the synthetically generated differentially expressed genes:

| Dataset | Method | Top-5 Accuracy | Top-10 Accuracy | Top-20 Accuracy |
|---------|--------|---------------|----------------|-----------------|
| Small (5 diff genes) | Quantum | 100% | N/A | N/A |
| Small (5 diff genes) | Classical | 100% | N/A | N/A |
| Medium (10 diff genes) | Quantum | 100% | 100% | N/A |
| Medium (10 diff genes) | Classical | 100% | 100% | N/A |
| Large (30 diff genes) | Quantum | 100% | 100% | 100% |
| Large (30 diff genes) | Classical | 100% | 100% | 100% |

Both methods achieved perfect accuracy in identifying the true differentially expressed genes across all dataset sizes. This demonstrates the robustness of both approaches for identifying genes with varying expression patterns along pseudotime.

## Gene Ranking and Prioritization

While both methods identified the same differentially expressed genes, they ranked them differently:

### Small Dataset (Top 5 genes)

**Quantum Method Top 5:**
1. DIFF_GENE_2 (Score: 4.18)
2. DIFF_GENE_4 (Score: 1.09)
3. DIFF_GENE_1 (Score: 1.02)
4. DIFF_GENE_3 (Score: 0.98)
5. DIFF_GENE_0 (Score: 0.73)

**Classical Method Top 5:**
1. DIFF_GENE_4 (Score: 0.94)
2. DIFF_GENE_1 (Score: 0.92)
3. DIFF_GENE_0 (Score: 0.92)
4. DIFF_GENE_3 (Score: 0.91)
5. DIFF_GENE_2 (Score: 0.52)

### Medium Dataset (Top 5 genes)

**Quantum Method Top 5:**
1. DIFF_GENE_5 (Score: 2.48)
2. DIFF_GENE_8 (Score: 2.30)
3. DIFF_GENE_2 (Score: 2.08)
4. DIFF_GENE_4 (Score: 0.64)
5. DIFF_GENE_7 (Score: 0.64)

**Classical Method Top 5:**
1. DIFF_GENE_1 (Score: 0.93)
2. DIFF_GENE_4 (Score: 0.93)
3. DIFF_GENE_7 (Score: 0.92)
4. DIFF_GENE_6 (Score: 0.90)
5. DIFF_GENE_0 (Score: 0.89)

## Overlap Analysis

We measured the overlap in top differentially expressed genes between methods:

| Dataset | Top-5 Overlap | Top-10 Overlap | Top-20 Overlap |
|---------|---------------|----------------|----------------|
| Small | 100% | N/A | N/A |
| Medium | 40% | 100% | N/A |
| Large | 0% | 0% | 50% |

These results reveal an interesting pattern: while both methods achieve perfect accuracy in identifying differentially expressed genes, they prioritize them differently. This divergence increases with dataset size:

- In the small dataset, both methods select the same top 5 genes (100% overlap)
- In the medium dataset, only 40% overlap in the top 5, but 100% overlap in the top 10
- In the large dataset, completely different genes in the top 5 and top 10, with only 50% overlap in the top 20

This suggests that each method is sensitive to different aspects of gene expression dynamics, potentially making them complementary approaches for gene discovery.

## Memory Usage

| Dataset | Memory Usage |
|---------|-------------|
| Small | 56.1MB |
| Medium | 62.8MB |
| Large | 78.4MB |

Memory usage remains quite efficient across all dataset sizes, with only a modest increase as the dataset size grows 25-fold (from 20 to 500 genes). This demonstrates the excellent scalability of the implementation.

## Pattern Detection Analysis

Analysis of gene expression patterns shows that both methods can detect various patterns:

1. **Linear patterns**: Both methods excel at identifying genes with linear expression changes
2. **Bell-curve patterns**: The quantum method shows particular strength in identifying genes with non-monotonic expression patterns, as seen in the top genes visualization
3. **Sigmoidal patterns**: Both methods can detect these transition patterns, with the quantum method often giving them higher ranking

The eigenvalue distribution from the quantum analysis (as seen in the `eigenvalues.png` visualization) shows a clear separation between the significant components driving gene expression changes and background variation, contributing to the method's high accuracy.

## How to Reproduce These Benchmarks

1. Clone the repository and install dependencies
```bash
git clone https://github.com/YourUsername/EntangleDE.git
cd EntangleDE
pip install -r requirements.txt
```

2. Run the benchmark script
```bash
# For small dataset
python benchmark.py --size small

# For medium dataset
python benchmark.py --size medium

# For large dataset
python benchmark.py --size large
```

## Conclusion

The EntangleDE quantum differential expression analysis method shows significant advantages over classical approaches:

1. **Dramatic scalability**: Execution time scales much better with the quantum approach, with speedups increasing from 6x (small dataset) to 42x (large dataset)

2. **Perfect accuracy**: Both methods achieve 100% accuracy in identifying differentially expressed genes, validating the quantum approach's effectiveness

3. **Complementary gene prioritization**: The quantum method ranks genes differently than classical approaches, particularly in larger datasets, suggesting it may detect different aspects of gene expression dynamics

4. **Efficient memory usage**: The quantum approach maintains modest memory requirements even as dataset size increases 25-fold

5. **Pattern recognition strengths**: The quantum method appears particularly adept at capturing non-monotonic and complex expression patterns

These results demonstrate that the quantum approach not only matches the accuracy of classical methods but dramatically outperforms them in computational efficiency. The different gene prioritization suggests that the methods are sensitive to different aspects of gene expression dynamics, making them potentially complementary techniques for gene discovery.

These benchmarks represent performance on simulation of quantum circuits using Qiskit. As quantum computing hardware improves, we expect the performance characteristics to improve even further.

## Future Work

1. Test on larger, real-world scRNA-seq datasets with varied biological conditions
2. Optimize the Hamiltonian embedding approach for even better performance
3. Explore different classical benchmarking methods 
4. Evaluate the methods on more complex expression patterns
5. Test with different parameter configurations to determine optimal settings for various dataset types
6. Investigate the biological significance of genes ranked differently by the two methods