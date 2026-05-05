# HLA-B39 Simulation Analysis

Analysis pipeline for molecular dynamics simulation data of HLA-B39 subtypes and B38:01 in peptide-loaded and peptide-free states.

## Overview

This repository analyzes protein interaction energies and network centrality changes in HLA class I molecules using molecular dynamics simulation data. The analysis focuses on identifying residues and interaction pairs that respond to peptide binding and differ between allelic variants.

## Repository Structure

### Analysis Scripts

#### [1_preprocessing.py](1_preprocessing.py)
Preprocessing pipeline that processes molecular dynamics simulation data from gRINN outputs.

**Steps:**
1. Load interaction energy data from gRINN outputs
2. Annotate residue pairs with structural information
3. Identify interaction pairs with consistent energies across replicates (ANOVA)
4. Compare peptide-loaded vs peptide-free systems (t-test)
5. Compare B*39 subtypes with B*38:01 reference
6. Construct protein energy networks (Ribeiro-Ortiz method)
7. Calculate betweenness centrality for all residues
8. Identify residues with significant centrality changes upon peptide binding

**Output:** 8 CSV files containing preprocessed datasets

**Runtime:** Several minutes (network construction is computationally intensive)

#### [2_analysis.ipynb](2_analysis.ipynb)
Jupyter notebook for visualization and exploratory analysis.

**Sections:**
- Summary statistics of interaction changes by allele and peptide
- 3D structure visualization of significantly affected residue pairs
- Conservation score integration with network centrality analysis
- Polymorphic position analysis across alleles

**Prerequisites:** Run `1_preprocessing.py` first to generate required CSV files

## Workflow

```bash
# Step 1: Run preprocessing (once)
python3 1_preprocessing.py

# Step 2: Run analysis (interactive)
jupyter notebook 2_analysis.ipynb
```

## Data Requirements

The preprocessing script expects gRINN outputs at:
```
BASE_FOLDER/*/grinn_output_skip10/energies_intEnVdW.csv
BASE_FOLDER/*/grinn_output_skip10/system_dry.pdb
```

## Methods

**Statistical Tests:**
- ANOVA: Identify consistent interactions/residues across replicates (p > 0.05)
- t-test: Compare peptide-loaded vs peptide-free conditions (p < 0.05)
- Energy threshold: Absolute mean difference > 2 kJ/mol

**Network Analysis:**
- Protein energy networks constructed using Ribeiro-Ortiz method
- Betweenness centrality calculated for all residues across simulation frames
- Equilibration: Frames > 20 retained for analysis

## Output Files

| File | Description |
|------|-------------|
| `intEnVdW_2025_07_10.csv` | Raw interaction energies with annotations |
| `df_cons_pairs_saved_2025_07_10.csv` | Consistent interaction pairs (ANOVA p > 0.05) |
| `df_sig_aff_pairs_loaded_2025_07_10.csv` | Pairs affected by peptide loading |
| `df_pairs_only_in_pmhc_2025_07_10.csv` | Pairs only present in peptide-loaded systems |
| `df_sig_aff_pairs_b39_2025_07_10.csv` | Pairs differing between B*39 and B*38:01 |
| `df_bc_equil_2025_07_10.csv` | Betweenness centrality values (equilibrated frames) |
| `df_bc_equil_cons_resids_2025_07_10.csv` | Residues with consistent BC across replicates |
| `df_bc_sigaff_resids_2026_01_27.csv` | Residues with significant BC changes |

## Dependencies

```
prody
pandas
numpy
scipy
matplotlib
seaborn
networkx
py3Dmol
natsort
tqdm
```
