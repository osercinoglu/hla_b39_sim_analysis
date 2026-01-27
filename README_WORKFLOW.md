# HLA-B39 Analysis Workflow

The analysis pipeline has been split into two separate files for easier management:

## Files Created

### 1. `1_preprocessing.py` (Python Script)
**Purpose:** Performs all computationally intensive preprocessing steps

**Sections:**
- 1.1: Load and Process Interaction Energies
- 1.2: Add Residue Annotations
- 1.3: Identify Consistent Interaction Pairs
- 1.4: Compare Peptide-Loaded vs Peptide-Free Systems
- 1.5: Compare B*39:01/06 Subtypes with B*38:01
- 1.6: Construct Protein Energy Networks
- 1.7: Identify Residues with Consistent Betweenness Centrality
- 1.8: Identify Residues with Significant BC Changes

**Output:** 8 CSV files containing preprocessed data

**Usage:**
```bash
python3 1_preprocessing.py
```

**Note:** This script takes several minutes to complete due to network construction and betweenness centrality calculations.

---

### 2. `2_analysis.ipynb` (Jupyter Notebook)
**Purpose:** Analysis and visualization using preprocessed data

**Sections:**
- 2.2: Summary Statistics (interaction changes by allele/peptide)
- 2.3: 3D Structure Visualization (attractive/repulsive interactions)
- 2.4: Conservation and Network Centrality Analysis
- 2.5: Polymorphic Position Analysis

**Prerequisites:** Must run `1_preprocessing.py` first

**Usage:**
```bash
jupyter notebook 2_analysis.ipynb
```

---

## Workflow

1. **Run preprocessing** (once):
   ```bash
   python3 1_preprocessing.py
   ```
   This generates all required CSV files.

2. **Run analysis** (as many times as needed):
   ```bash
   jupyter notebook 2_analysis.ipynb
   ```
   The notebook loads data from CSV files - no need to re-run preprocessing.

---

## Generated CSV Files

The preprocessing script generates the following files:

1. `intEnVdW_2025_07_10.csv` - Interaction energy data
2. `df_cons_pairs_saved_2025_07_10.csv` - Consistent interaction pairs
3. `df_sig_aff_pairs_loaded_2025_07_10.csv` - Significantly affected pairs (peptide loading)
4. `df_pairs_only_in_pmhc_2025_07_10.csv` - Pairs only in peptide-loaded systems
5. `df_sig_aff_pairs_b39_2025_07_10.csv` - B*39 vs B*38:01 comparisons
6. `df_bc_equil_2025_07_10.csv` - Betweenness centrality equilibrium data
7. `df_bc_equil_cons_resids_2025_07_10.csv` - Residues with consistent BC
8. `df_bc_sigaff_resids_2026_01_27.csv` - Residues with significant BC changes

All files are loaded automatically by the analysis notebook.

---

## Benefits of This Structure

✓ **Separation of concerns:** Preprocessing and analysis are clearly separated
✓ **Time savings:** No need to re-run lengthy preprocessing for analysis changes
✓ **Easier debugging:** Python script can be run from command line with print statements
✓ **Better organization:** Clear execution flow with no cell ordering confusion
✓ **Reproducibility:** Preprocessing runs the same way every time

