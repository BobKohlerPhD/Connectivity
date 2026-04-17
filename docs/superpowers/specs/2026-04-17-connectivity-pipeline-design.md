# Connectivity Pipeline Design Specification (2026-04-17)

## Overview
A configuration-driven, modular R pipeline for Connectome-based Predictive Modeling (CPM). This system transforms a collection of fragmented scripts into a robust, "state-of-the-field" toolkit for neuroimaging connectivity analysis.

## Goals
- **Modularity:** Decoupled functional units for loading, cleaning, modeling, and plotting.
- **Reproducibility:** All parameters and paths managed via a central `config.yaml`.
- **Atlas-Agnosticism:** Abstracted mapping logic for easy integration of any brain atlas (Shen, Glasser, Schaefer, etc.).
- **Robustness:** Strict data validation (symmetry, NaN checks) and error logging.
- **User-Friendliness:** Single CLI entry point (`main.R`) for end-to-end execution.

## Architecture

### Directory Structure
```text
Connectivity/
├── config/
│   ├── default_config.yaml      # Master configuration
│   └── atlas_shen.yaml          # Shen-specific mapping & labels
├── R/                           # Modularized core logic
│   ├── io_utils.R               # Loading/Saving, format conversion
│   ├── preprocessing.R          # Cleaning, alignment, covariate handling
│   ├── cpm_engine.R             # Core CPM modeling and CV logic
│   ├── posthoc_analysis.R       # Overlap stats, external mask testing
│   └── plotting_utils.R         # Unified heatmap and network plotting
├── main.R                       # CLI entry point
├── atlases/                     # Atlas configuration and mapping templates
│   ├── Shen/                    # Original Shen files (NodeList, canonical nodes, etc.)
│   ├── Glasser/                 # Glasser (HCP) mapping (360 nodes)
│   ├── Schaefer/                # Schaefer mapping (e.g., 200, 400 node versions)
│   └── Power/                   # Power mapping (264 nodes)
└── docs/                        # Technical documentation

## Common Atlas Support
The pipeline includes pre-configured templates for common research atlases. Each atlas requires a standardized `mapping.yaml` file that defines:
- **Node Order:** For canonical network rearrangement.
- **Network Labels:** For visualization and summary reports.
- **Node-to-Network Mapping:** To aggregate edges into functional modules.

### Supported Atlas Templates
- **Shen (268 nodes):** Default mapping provided.
- **Glasser (360 nodes):** Multi-modal parcellation (HCP).
- **Schaefer (User-defined nodes):** Support for 100, 200, 400, etc., node versions.
- **Power (264 nodes):** Common functional parcellation.
- **Gordon (333 nodes):** Resting-state functional parcellation.

## Atlas Mapping Schema (`atlases/<atlas_name>_mapping.yaml`)
```yaml
name: "Glasser"
num_nodes: 360
network_labels:
  - "Primary Visual"
  - "Dorsal Stream"
  - "Frontoparietal"
  - "Default Mode"
  # ... rest of labels
node_order: [1, 5, 20, 35, ...] # Indices for rearrangement
network_mapping: "atlases/Glasser/node_to_network.csv" # CSV mapping individual nodes
```
```

### Core Modules

#### 1. `io_utils.R`
- **Standardized Loading:** Supports `.mat`, `.rds`, and `.csv` via `R.matlab`, `readr`, and `tidyverse`.
- **Validation:** Automatic checks for matrix symmetry, node count consistency, and subject ID alignment.
- **Logging:** Descriptive errors for missing data or malformed files.

#### 2. `preprocessing.R`
- **Alignment:** Computes the "perfect intersection" between connectivity matrices and behavioral metadata.
- **Cleaning:** Removes subjects with `NaN`, `Inf`, or missing edges.
- **Covariates:** Standardized dummy coding and partial correlation setup for regression.

#### 3. `cpm_engine.R`
- **Modeling:** Implements Linear, Ridge, and Logistic regression modes.
- **Cross-Validation:** Robust k-fold implementation with parallel execution using the `future` package.
- **Consensus:** Aggregates predictive edges into consensus networks with configurable selection thresholds.

#### 4. `posthoc_analysis.R`
- **Overlap Stats:** Hypergeometric, Jaccard, and Cosine similarity for two- and three-way overlaps.
- **External Mask Testing:** Evaluates discovered networks on new data using `lm`, `glm`, or `lmer`.

#### 5. `plotting_utils.R`
- **Heatmaps:** Publication-ready `ggplot2` visualizations.
- **Metadata-Driven:** Decouples plotting from atlas structure; uses external mapping files for node order and labels.

## Configuration System (`config.yaml`)
Users control the pipeline via a YAML file, avoiding the need to edit core R logic.

```yaml
project_name: "Connectivity_Study_2024"
output_dir: "./results"

data:
  matrices: "data/raw/connectivity.rds"
  outcomes: "data/raw/behavioral_data.csv"
  subject_id_var: "sub_id"
  outcome_var: "total_score"
  covariates: ["age", "sex", "motion_params"]

cpm_params:
  k_folds: 10
  iterations: 100
  p_threshold: 0.05
  mode: "linear" # options: linear, ridge, logistic
  z_score: true

atlas:
  name: "Shen"
  mapping_file: "atlases/shen_mapping.yaml"
```

## CLI Usage
```bash
Rscript main.R --config config/my_study.yaml
```

## Error Handling & Robustness
- **Logging:** All pipeline steps log to both the terminal and a `run.log` file.
- **Resiliency:** Intermediate results are saved to allow for recovery in long-running jobs.
- **Validation:** Halts immediately if critical data integrity checks fail (e.g., non-symmetric matrices).

## Testing Plan
- **Unit Tests:** Validate individual functions in `io_utils`, `preprocessing`, and `cpm_engine`.
- **Integration Tests:** End-to-end runs with synthetic datasets to verify pipeline flow.
- **Reproducibility:** Seeded runs to ensure bit-perfect result matching.
