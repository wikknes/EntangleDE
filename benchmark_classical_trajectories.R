#!/usr/bin/env Rscript
# benchmark_classical_trajectories.R
# Script to benchmark classical trajectory inference methods: Monocle3 and Slingshot
# Usage: Rscript benchmark_classical_trajectories.R --dataset small

# Function to check and install packages
check_and_install <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(paste0("Package ", pkg_name, " not installed. Attempting to install...\n"))
    install.packages(pkg_name, repos = "https://cloud.r-project.org/")
    
    # Check if BioConductor package
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      if (pkg_name %in% c("SingleCellExperiment", "slingshot")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org/")
        }
        BiocManager::install(pkg_name)
      }
    }
    
    # Check if still not installed
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      stop(paste0("Failed to install ", pkg_name, ". Please install manually."))
    }
  }
  cat(paste0(pkg_name, " is installed.\n"))
}

# Check and install required packages
required_packages <- c("optparse", "SingleCellExperiment", "reticulate", 
                      "ggplot2", "gridExtra", "reshape2", "plyr", 
                      "dplyr", "readr", "Seurat", "slingshot")

# Monocle3 requires special handling
if (!requireNamespace("monocle3", quietly = TRUE)) {
  cat("monocle3 not installed. It requires manual installation.\n")
  cat("Please install using instructions at: https://cole-trapnell-lab.github.io/monocle3/docs/installation/\n")
  cat("Continuing without monocle3 benchmarks...\n")
  run_monocle <- FALSE
} else {
  run_monocle <- TRUE
  cat("monocle3 is installed.\n")
}

for (pkg in required_packages) {
  check_and_install(pkg)
}

# Load packages
library(optparse)
if (run_monocle) library(monocle3)
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)
library(dplyr)
library(readr)

# Setup command line options
option_list <- list(
  make_option(c("--dataset"), type="character", default="small",
              help="Dataset size to use [small, medium, large]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Load Python environment to use the same synthetic data
np <- import("numpy")
pd <- import("pandas")

# Set paths based on dataset size
size <- opt$dataset
data_dir <- file.path("benchmark_trajectory_results", size)
output_dir <- file.path("benchmark_classical_results", size)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load synthetic data
read_synthetic_data <- function(size) {
  # Read expression data
  expr_file <- file.path(data_dir, "synthetic_expression.csv")
  expr_data <- read.csv(expr_file, row.names = 1)
  
  # Read pseudotime data
  pseudotime_file <- file.path(data_dir, "synthetic_pseudotime.csv")
  pseudotime_data <- read.csv(pseudotime_file)
  
  # Read branch data
  branch_file <- file.path(data_dir, "synthetic_branches.csv")
  branch_data <- read.csv(branch_file)
  
  return(list(
    expression = as.matrix(expr_data),
    pseudotime = pseudotime_data$pseudotime,
    branches = branch_data$branch
  ))
}

# Load data
print(paste("Loading", size, "synthetic dataset..."))
data <- read_synthetic_data(size)

# Convert to cell x gene format for R tools
expr_matrix <- t(data$expression)
rownames(expr_matrix) <- paste0("Cell_", 1:nrow(expr_matrix))
colnames(expr_matrix) <- paste0("Gene_", 1:ncol(expr_matrix))

# Create metadata
cell_metadata <- data.frame(
  cell_id = rownames(expr_matrix),
  true_pseudotime = data$pseudotime,
  true_branch = data$branches,
  row.names = rownames(expr_matrix)
)

# Print dataset info
print(paste("Dataset size:", size))
print(paste("Cells:", nrow(expr_matrix)))
print(paste("Genes:", ncol(expr_matrix)))
print(paste("Branches:", length(unique(data$branches))))

# Helper function to compute Kendall's Tau
compute_kendall <- function(true_time, inferred_time) {
  result <- cor.test(true_time, inferred_time, method = "kendall")
  return(list(
    tau = result$estimate,
    p.value = result$p.value
  ))
}

# Helper function to compute Adjusted Rand Index
compute_ari <- function(true_clusters, inferred_clusters) {
  require(mclust)
  return(adjustedRandIndex(true_clusters, inferred_clusters))
}

# Helper function for silhouette score
compute_silhouette <- function(data, clusters) {
  require(cluster)
  if (length(unique(clusters)) <= 1) {
    return(0)
  }
  sil <- silhouette(as.numeric(clusters), dist(data))
  return(mean(sil[,3]))
}

# Run Monocle3
run_monocle3 <- function(expr_matrix, cell_metadata) {
  print("Running Monocle3...")
  start_time <- Sys.time()
  
  # Create CDS object
  gene_metadata <- data.frame(
    gene_id = colnames(expr_matrix),
    gene_short_name = colnames(expr_matrix),
    row.names = colnames(expr_matrix)
  )
  
  cds <- new_cell_data_set(
    expression_data = as.matrix(expr_matrix),
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  
  # Preprocess
  cds <- preprocess_cds(cds, num_dim = 20)
  
  # Reduce dimension
  cds <- reduce_dimension(cds)
  
  # Cluster cells
  cds <- cluster_cells(cds)
  
  # Learn graph
  cds <- learn_graph(cds)
  
  # Order cells
  # Use the cell with the lowest true pseudotime as the root
  root_cell <- rownames(cell_metadata)[which.min(cell_metadata$true_pseudotime)]
  cds <- order_cells(cds, root_cells = root_cell)
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract pseudotime and clusters
  pseudotime <- pseudotime(cds)
  clusters <- clusters(cds)
  
  # Compute metrics
  kendall_result <- compute_kendall(cell_metadata$true_pseudotime, pseudotime)
  ari <- compute_ari(cell_metadata$true_branch, as.numeric(clusters))
  sil <- compute_silhouette(reducedDim(cds, "UMAP"), as.numeric(clusters))
  
  # Plot results
  p1 <- ggplot(data.frame(
    UMAP1 = reducedDim(cds, "UMAP")[,1],
    UMAP2 = reducedDim(cds, "UMAP")[,2],
    Pseudotime = pseudotime,
    Cluster = clusters
  )) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = Pseudotime)) +
    scale_color_viridis_c() +
    ggtitle("Monocle3 Pseudotime") +
    theme_minimal()
  
  p2 <- ggplot(data.frame(
    UMAP1 = reducedDim(cds, "UMAP")[,1],
    UMAP2 = reducedDim(cds, "UMAP")[,2],
    Pseudotime = pseudotime,
    Cluster = clusters
  )) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    ggtitle("Monocle3 Clusters") +
    theme_minimal()
  
  # Save plots
  g <- grid.arrange(p1, p2, ncol = 2)
  ggsave(file.path(output_dir, "monocle3_results.png"), g, width = 10, height = 5)
  
  return(list(
    execution_time = execution_time,
    pseudotime = pseudotime,
    clusters = clusters,
    kendall_tau = kendall_result$tau,
    p_value = kendall_result$p.value,
    ari = ari,
    silhouette = sil,
    cds = cds
  ))
}

# Run Slingshot
run_slingshot <- function(expr_matrix, cell_metadata) {
  print("Running Slingshot...")
  start_time <- Sys.time()
  
  # Create Seurat object for preprocessing
  seu <- CreateSeuratObject(counts = t(expr_matrix))
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 20)
  seu <- RunUMAP(seu, dims = 1:20)
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.8)
  
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seu)
  
  # Use the cell with lowest true pseudotime as starting cluster
  root_cell <- which.min(cell_metadata$true_pseudotime)
  root_cluster <- seu$seurat_clusters[root_cell]
  
  # Run Slingshot
  sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', start.clus = root_cluster)
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Extract pseudotime and clusters from the first lineage (simplification)
  if ("slingPseudotime_1" %in% colnames(colData(sce))) {
    pseudotime <- sce$slingPseudotime_1
  } else {
    # If only one lineage
    pseudotime <- sce$slingPseudotime
  }
  
  # Handle NA values in pseudotime by assigning them max value
  # (cells not assigned to this lineage are likely later in development)
  pseudotime[is.na(pseudotime)] <- max(pseudotime, na.rm = TRUE)
  
  clusters <- seu$seurat_clusters
  
  # Compute metrics
  kendall_result <- compute_kendall(cell_metadata$true_pseudotime, pseudotime)
  ari <- compute_ari(cell_metadata$true_branch, as.numeric(clusters))
  sil <- compute_silhouette(reducedDims(sce)$UMAP, as.numeric(clusters))
  
  # Plot results
  sce$pseudotime <- pseudotime
  sce$clusters <- clusters
  
  df <- data.frame(
    UMAP1 = reducedDims(sce)$UMAP[,1],
    UMAP2 = reducedDims(sce)$UMAP[,2],
    Pseudotime = pseudotime,
    Cluster = clusters
  )
  
  p1 <- ggplot(df) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = Pseudotime)) +
    scale_color_viridis_c() +
    ggtitle("Slingshot Pseudotime") +
    theme_minimal()
  
  p2 <- ggplot(df) +
    geom_point(aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    ggtitle("Slingshot Clusters") +
    theme_minimal()
  
  # Save plots
  g <- grid.arrange(p1, p2, ncol = 2)
  ggsave(file.path(output_dir, "slingshot_results.png"), g, width = 10, height = 5)
  
  return(list(
    execution_time = execution_time,
    pseudotime = pseudotime,
    clusters = clusters,
    kendall_tau = kendall_result$tau,
    p_value = kendall_result$p.value,
    ari = ari,
    silhouette = sil,
    sce = sce
  ))
}

# Run benchmarks
results_list <- list()

# Run Monocle3 if available
if (run_monocle) {
  monocle3_results <- try(run_monocle3(expr_matrix, cell_metadata))
  if (!inherits(monocle3_results, "try-error")) {
    results_list[["monocle3"]] <- monocle3_results
  } else {
    cat("Error running Monocle3 benchmark. Skipping...\n")
  }
} else {
  cat("Skipping Monocle3 benchmarks as it's not installed.\n")
}

# Run Slingshot
slingshot_results <- try(run_slingshot(expr_matrix, cell_metadata))
if (!inherits(slingshot_results, "try-error")) {
  results_list[["slingshot"]] <- slingshot_results
} else {
  cat("Error running Slingshot benchmark. Skipping...\n")
}

# Summarize results
print("===== BENCHMARK RESULTS =====")
print(paste("Dataset:", size))

# Function to safely access results
get_metric <- function(results, metric, default = NA) {
  if (is.null(results) || inherits(results, "try-error") || is.null(results[[metric]])) {
    return(default)
  }
  return(results[[metric]])
}

print("\nExecution Time:")
if ("monocle3" %in% names(results_list)) {
  print(paste("Monocle3:", round(get_metric(results_list$monocle3, "execution_time"), 2), "seconds"))
} else {
  print("Monocle3: Not run")
}
if ("slingshot" %in% names(results_list)) {
  print(paste("Slingshot:", round(get_metric(results_list$slingshot, "execution_time"), 2), "seconds"))
} else {
  print("Slingshot: Not run")
}

print("\nPseudotime Correlation (Kendall's Tau):")
if ("monocle3" %in% names(results_list)) {
  print(paste("Monocle3:", round(get_metric(results_list$monocle3, "kendall_tau"), 3), 
              "p-value:", round(get_metric(results_list$monocle3, "p_value"), 3)))
} else {
  print("Monocle3: Not run")
}
if ("slingshot" %in% names(results_list)) {
  print(paste("Slingshot:", round(get_metric(results_list$slingshot, "kendall_tau"), 3), 
              "p-value:", round(get_metric(results_list$slingshot, "p_value"), 3)))
} else {
  print("Slingshot: Not run")
}

print("\nClustering Accuracy (ARI):")
if ("monocle3" %in% names(results_list)) {
  print(paste("Monocle3:", round(get_metric(results_list$monocle3, "ari"), 3)))
} else {
  print("Monocle3: Not run")
}
if ("slingshot" %in% names(results_list)) {
  print(paste("Slingshot:", round(get_metric(results_list$slingshot, "ari"), 3)))
} else {
  print("Slingshot: Not run")
}

print("\nCluster Coherence (Silhouette):")
if ("monocle3" %in% names(results_list)) {
  print(paste("Monocle3:", round(get_metric(results_list$monocle3, "silhouette"), 3)))
} else {
  print("Monocle3: Not run")
}
if ("slingshot" %in% names(results_list)) {
  print(paste("Slingshot:", round(get_metric(results_list$slingshot, "silhouette"), 3)))
} else {
  print("Slingshot: Not run")
}

# Prepare results dataframe
tool_names <- character()
execution_times <- numeric()
kendall_taus <- numeric()
p_values <- numeric()
aris <- numeric()
silhouettes <- numeric()

# Add results for each available method
if ("monocle3" %in% names(results_list)) {
  tool_names <- c(tool_names, "Monocle3")
  execution_times <- c(execution_times, get_metric(results_list$monocle3, "execution_time"))
  kendall_taus <- c(kendall_taus, get_metric(results_list$monocle3, "kendall_tau"))
  p_values <- c(p_values, get_metric(results_list$monocle3, "p_value"))
  aris <- c(aris, get_metric(results_list$monocle3, "ari"))
  silhouettes <- c(silhouettes, get_metric(results_list$monocle3, "silhouette"))
}

if ("slingshot" %in% names(results_list)) {
  tool_names <- c(tool_names, "Slingshot")
  execution_times <- c(execution_times, get_metric(results_list$slingshot, "execution_time"))
  kendall_taus <- c(kendall_taus, get_metric(results_list$slingshot, "kendall_tau"))
  p_values <- c(p_values, get_metric(results_list$slingshot, "p_value"))
  aris <- c(aris, get_metric(results_list$slingshot, "ari"))
  silhouettes <- c(silhouettes, get_metric(results_list$slingshot, "silhouette"))
}

# Create dataframe if we have results
if (length(tool_names) > 0) {
  results_df <- data.frame(
    Tool = tool_names,
    Execution_Time_Seconds = execution_times,
    Kendall_Tau = kendall_taus,
    P_Value = p_values,
    ARI = aris,
    Silhouette = silhouettes
  )
} else {
  results_df <- data.frame(
    Tool = character(),
    Execution_Time_Seconds = numeric(),
    Kendall_Tau = numeric(),
    P_Value = numeric(),
    ARI = numeric(),
    Silhouette = numeric()
  )
  cat("No results available to save.\n")
}

write.csv(results_df, file.path(output_dir, paste0(size, "_classical_metrics.csv")), row.names = FALSE)

# Compare with EntangleDE and Scanpy
try({
  q_metrics_file <- file.path("benchmark_trajectory_results", paste0(size, "_trajectory_metrics.csv"))
  if (file.exists(q_metrics_file)) {
    q_metrics <- read.csv(q_metrics_file)
    
    # Create base comparison dataframe with EntangleDE and Scanpy
    tool_names <- c("EntangleDE", "Scanpy")
    execution_times <- c(q_metrics$quantum_execution_time, q_metrics$scanpy_execution_time)
    kendall_taus <- c(q_metrics$quantum_kendall_tau, q_metrics$scanpy_kendall_tau)
    aris <- c(q_metrics$quantum_ari, q_metrics$scanpy_ari)
    silhouettes <- c(q_metrics$quantum_silhouette, q_metrics$scanpy_silhouette)
    
    # Add Monocle3 if available
    if ("monocle3" %in% names(results_list)) {
      tool_names <- c(tool_names, "Monocle3")
      execution_times <- c(execution_times, get_metric(results_list$monocle3, "execution_time"))
      kendall_taus <- c(kendall_taus, get_metric(results_list$monocle3, "kendall_tau"))
      aris <- c(aris, get_metric(results_list$monocle3, "ari"))
      silhouettes <- c(silhouettes, get_metric(results_list$monocle3, "silhouette"))
    }
    
    # Add Slingshot if available
    if ("slingshot" %in% names(results_list)) {
      tool_names <- c(tool_names, "Slingshot")
      execution_times <- c(execution_times, get_metric(results_list$slingshot, "execution_time"))
      kendall_taus <- c(kendall_taus, get_metric(results_list$slingshot, "kendall_tau"))
      aris <- c(aris, get_metric(results_list$slingshot, "ari"))
      silhouettes <- c(silhouettes, get_metric(results_list$slingshot, "silhouette"))
    }
    
    comparison_df <- data.frame(
      Tool = tool_names,
      Execution_Time = execution_times,
      Kendall_Tau = kendall_taus,
      ARI = aris,
      Silhouette = silhouettes
    )
  
  # Save comprehensive comparison
  write.csv(comparison_df, file.path(output_dir, paste0(size, "_all_methods_comparison.csv")), row.names = FALSE)
  
  # Create comparison plots
  # Execution time comparison
  p1 <- ggplot(comparison_df, aes(x = Tool, y = Execution_Time, fill = Tool)) +
    geom_bar(stat = "identity") +
    labs(title = "Execution Time Comparison", y = "Time (seconds)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Pseudotime accuracy
  p2 <- ggplot(comparison_df, aes(x = Tool, y = Kendall_Tau, fill = Tool)) +
    geom_bar(stat = "identity") +
    labs(title = "Pseudotime Accuracy (Kendall's Tau)", y = "Kendall's Tau") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Clustering accuracy
  p3 <- ggplot(comparison_df, aes(x = Tool, y = ARI, fill = Tool)) +
    geom_bar(stat = "identity") +
    labs(title = "Clustering Accuracy (ARI)", y = "Adjusted Rand Index") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plots
  g <- grid.arrange(p1, p2, p3, ncol = 2)
  ggsave(file.path(output_dir, "all_methods_comparison.png"), g, width = 12, height = 8)
  
  print("Created comprehensive comparison with EntangleDE and Scanpy")
})

print(paste("Results saved to", output_dir))
print("Benchmark complete!")