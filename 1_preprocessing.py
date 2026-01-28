"""
HLA-B39 Simulation Analysis - Preprocessing Script

This script performs all preprocessing steps:
- Load and process interaction energies from gRINN outputs
- Identify consistent interaction pairs across replicates
- Compare peptide-loaded vs peptide-free systems
- Compare B*39 subtypes with B*38:01
- Construct protein energy networks
- Calculate betweenness centrality
- Identify residues with significant changes

Output: Multiple CSV files for downstream analysis
"""

# Imports
import glob
from prody import *
import pandas as pd
import re
import numpy as np
import tqdm
import natsort
import seaborn as sns
import scipy.stats
import matplotlib.pyplot as plt
import py3Dmol
import networkx as nx

# Limits of p values:
# We do that because some p values can be reported to be zero.
np.finfo('float')


print("="*80)
print("HLA-B39 PREPROCESSING PIPELINE")
print("="*80)

print("\n[1/8] Loading gRINN outputs and processing interaction energies...")
BASE_FOLDER = '/mnt/MORIA/experiments/hla_b39/2022_10_14_bioe_leviathan_backup/T1DM_B39_traj_analysis/gmx_traj_data/'
grinn_outs = glob.glob(BASE_FOLDER+'*/grinn_output_skip10/energies_intEnVdW.csv')
dry_pdbs = glob.glob(BASE_FOLDER+'*/grinn_output_skip10/system_dry.pdb')

grinn_outs = [out for out in grinn_outs if 'B3801_ALL_2' not in out]

dry_pdbs = [pdb for pdb in dry_pdbs if 'B3801_ALL_2' not in pdb]

print(f"Found {len(grinn_outs)} gRINN output files")
print(f"Found {len(dry_pdbs)} dry PDB files")

zipped = zip(dry_pdbs, grinn_outs)
df_list = list()

# Define the theoretical frame columns
frame_cols = list(map(str,np.arange(0, 251, 1)))

print("Processing simulation systems...")
for zipd in tqdm.tqdm(list(zipped), desc="Loading systems"):
    dry_pdb = zipd[0]
    matches = re.search(r'.*gmx_traj_data/(B\d+)(?:_([^/]*?))?_(\d+)/grinn_output_skip10/',dry_pdb)
    if matches:
        allele = matches.groups()[0]
        peptide = matches.groups()[1] if matches.groups()[1] else ''
        replica = int(matches.groups()[2])

    df_intEn = pd.read_csv(zipd[1])
    
    # Check which frame columns actually exist in the CSV
    available_frame_cols = [col for col in frame_cols if col in df_intEn.columns]
    
    df_intEn['repeat'] = replica
    df_intEn['allele'] = allele
    df_intEn['peptide'] = peptide
    df_intEn['system'] = allele + '_' + peptide if peptide else allele
    df_intEn['dry_pdb_file'] = dry_pdb

    # Use only the available frame columns
    cols_to_select = ['Pair_indices','repeat','allele','peptide','dry_pdb_file','system'] + available_frame_cols
    df_intEn = df_intEn[cols_to_select]
    
    # Change column names from "Unnamed: 0" to "pair"
    df_intEn.rename({'Pair_indices': 'pair'}, axis=1, inplace=True)

    # Limit to 377 columns as originally intended, but only if we have that many
    max_cols = min(377, len(df_intEn.columns))
    df_intEn = df_intEn[df_intEn.columns[:max_cols]]
    
    df_list.append(df_intEn)

df_intEn = pd.concat(df_list)

print(f"Combined {len(df_list)} dataframes into one")
print(f"Total interaction pairs: {len(df_intEn)}")
print(f"Total columns before dropna: {len(df_intEn.columns)}")

df_intEn.dropna(axis=1, inplace=True)

print(f"Total columns after dropna: {len(df_intEn.columns)}")
print("Saving to intEnVdW_2025_07_10.csv...")

df_intEn.to_csv('intEnVdW_2025_07_10.csv')
print("✓ Saved intEnVdW_2025_07_10.csv")


################################################################################
# ## 1.2 Add Residue Annotations
#
# Add residue numbering and chain information to interaction pairs.
################################################################################

print("\n[2/8] Adding residue annotations to interaction pairs...")
df_intEn = pd.read_csv('intEnVdW_2025_07_10.csv')
print(f"Loaded {len(df_intEn)} interaction pairs")

dry_pdb_file = ''

for row in tqdm.tqdm(df_intEn.iterrows(), total=len(df_intEn), desc="Annotating residues"):
    pair = row[1]['pair']
    dry_pdb_file_row = row[1]['dry_pdb_file']
    if dry_pdb_file_row != dry_pdb_file:
        syst_dry = parsePDB(dry_pdb_file_row)
        dry_pdb_file = dry_pdb_file_row
    
    pairs = pair.split('-')
    pair1 = syst_dry.select(f'resindex {pairs[0]}')
    pair2 = syst_dry.select(f'resindex {pairs[1]}')
    resnum1 = pair1.getResnums()[0]
    resnum2 = pair2.getResnums()[0]
    chain1 = pair1.getChids()[0]
    chain2 = pair2.getChids()[0]
    df_intEn.loc[row[0],'resnum1'] = resnum1
    df_intEn.loc[row[0],'chain1'] = chain1
    df_intEn.loc[row[0],'resnum2'] = resnum2
    df_intEn.loc[row[0],'chain2'] = chain2
    resnum_chain1 = f'{resnum1}{chain1}'
    resnum_chain2 = f'{resnum2}{chain2}'
    resnum_chain12 = [resnum_chain1,resnum_chain2]
    resnum_chain12 = natsort.natsorted(resnum_chain12)
    resnum_chain12 = '-'.join(resnum_chain12)
    df_intEn.loc[row[0],'resnum_chain1'] = resnum_chain1
    df_intEn.loc[row[0],'resnum_chain2'] = resnum_chain2
    df_intEn.loc[row[0],'resnum_chain12'] = resnum_chain12


################################################################################
# ## 1.3 Identify Consistent Interaction Pairs
#
# Find residue pairs with stable interactions across simulation replicates (ANOVA p > 0.05).
################################################################################

print("\n[3/8] Identifying consistent interaction pairs across replicates...")
uniq_pairs = df_intEn['resnum_chain12'].unique()
print(f"Total unique pairs to analyze: {len(uniq_pairs)}")

consistent_pairs = list()
other_pairs = list()

p_values_cons = list()
p_values_other = list()

for pair in tqdm.tqdm(uniq_pairs, desc="Testing pair consistency"):
    df_intEn_pair = df_intEn[df_intEn['resnum_chain12']==pair].copy()

    resnum1 = df_intEn_pair['resnum1'].unique()[0]
    resnum2 = df_intEn_pair['resnum2'].unique()[0]
    chain1 = df_intEn_pair['chain1'].unique()[0]
    chain2 = df_intEn_pair['chain2'].unique()[0]

    if chain1 == chain2 and np.abs(resnum1-resnum2) < 4:
        # Don't take pairs that are too close along the covalently-bonded polypeptide chain into account.
        continue

    systems = df_intEn_pair['system'].unique()
    for system in systems:
        df_intEn_pair_system = df_intEn_pair[df_intEn_pair['system'] == system].copy()
        if not len(df_intEn_pair_system) in [3,4]:
            # Then we exclude this pair, it was not observed in at least two repeats 
            continue
        
        # Remove columns and prepare data
        columns_to_drop = ['Unnamed: 0','dry_pdb_file','system','pair','resnum1','chain1','resnum2','chain2',
                          'resnum_chain1','resnum_chain2','resnum_chain12','repeat','allele','peptide']
        
        # Only drop columns that actually exist
        columns_to_drop = [col for col in columns_to_drop if col in df_intEn_pair_system.columns]
        df_intEn_pair_system.drop(columns_to_drop, axis=1, inplace=True)
        
        # Check if we have any data left
        if df_intEn_pair_system.empty or df_intEn_pair_system.shape[1] == 0:
            continue
            
        df_values = df_intEn_pair_system.T
        
        # Check if we have enough replicates (rows) for comparison
        if len(df_values.columns) < 2:
            continue
            
        df_values.columns = np.arange(len(df_values.columns))

        # Calculate means for each replicate
        means = []
        for i in range(len(df_values.columns)):
            col_data = df_values[i].dropna()  # Remove NaN values
            if len(col_data) > 0:
                means.append(np.abs(np.mean(col_data)))
            else:
                means.append(0)
        
        means = np.array(means)
        
        if np.all(means > 0.8):
            try:
                # Prepare data for ANOVA
                data_for_anova = []
                for i in range(len(df_values.columns)):
                    col_data = df_values[i].dropna().values
                    if len(col_data) > 0:
                        data_for_anova.append(col_data)
                
                # Need at least 2 groups for ANOVA
                if len(data_for_anova) >= 2:
                    f_statistic, p_value = scipy.stats.f_oneway(*data_for_anova)
                    
                    if p_value > 0.05:
                        consistent_pairs.append(pair)
                        p_values_cons.append(p_value)
                    else:
                        other_pairs.append(pair)
                        p_values_other.append(p_value)
                        
            except Exception as e:
                print(f"Error with pair {pair} in system {system}: {e}")
                print(f"Data shapes: {[len(arr) for arr in data_for_anova] if 'data_for_anova' in locals() else 'No data'}")
                continue


df_cons_pairs = pd.DataFrame(consistent_pairs, columns=['resnum_chain12'])
df_cons_pairs['p_value'] = p_values_cons

df_other_pairs = pd.DataFrame(other_pairs, columns=['resnum_chain12'])
df_other_pairs['p_value'] = p_values_other

print(f"Found {len(consistent_pairs)} consistent pairs (p > 0.05)")
print(f"Found {len(other_pairs)} variable pairs (p < 0.05)")

df_cons_pairs = pd.merge(df_intEn, df_cons_pairs, how='inner', on='resnum_chain12')
df_other_pairs = pd.merge(df_intEn, df_other_pairs, how='inner', on='resnum_chain12')

# Define the final best list of frame columns
# These are the columms detected for each trajectory
frame_cols = list(map(str,np.arange(0, 239, 1)))
cols2_save = ['repeat','system','allele','peptide','resnum_chain12','p_value','resnum_chain1','resnum_chain2'] + frame_cols
print("Saving consistent pairs to df_cons_pairs_saved_2025_07_10.csv...")
df_cons_pairs[cols2_save].to_csv('df_cons_pairs_saved_2025_07_10.csv', index=False)
print("✓ Saved df_cons_pairs_saved_2025_07_10.csv")

pairs = df_cons_pairs['resnum_chain12'].unique()
print(f"\nAnalyzing {len(pairs)} consistent pairs for system differences...")

lonely_consistent_pairs = dict()

sig_aff_pairs = list()

for pair in tqdm.tqdm(pairs, desc="Testing system differences"):
    df_intEn_pair = df_cons_pairs[df_cons_pairs['resnum_chain12']==pair].copy()

    unique_systems = df_intEn_pair['system'].unique()

    if len(unique_systems) == 1:
        # Then this pair was consistently observed in only one simulation.
        # Add that to lonely_consistent_pairs, but don't do anything further in the loop.
        lonely_consistent_pairs[unique_systems[0]] = pair
        continue

    pairvals = dict()
    for system in unique_systems:
        df_pair_system = df_intEn_pair[df_intEn_pair['system']==system].copy()

        val_arrays = list()
        for row in df_pair_system.iterrows():
            val_arrays.append(np.array(row[1][frame_cols]))

        val_arrays = np.hstack(val_arrays)

        pairvals[system] = val_arrays

    f_statistic, p_value = scipy.stats.f_oneway(*[v for v in pairvals.values()])
    if p_value < 0.05:
        #print(f'Significantly affected pair: {pair} with p-value {p_value}')
        sig_aff_pairs.append(pair)

print(f"Found {len(sig_aff_pairs)} significantly affected pairs across systems")

sig_aff_pairs = [el.split('-') for el in sig_aff_pairs]

df_sig_aff_pairs = pd.DataFrame(sig_aff_pairs, columns=['rc1','rc2'])
df_sig_aff_res = np.hstack([df_sig_aff_pairs['rc1'].values,df_sig_aff_pairs['rc2'].values])
df_sig_aff_res = pd.DataFrame(df_sig_aff_res, columns=['res'])
df_sig_aff_res.value_counts()[:50]


################################################################################
# ## 1.4 Compare Peptide-Loaded vs Peptide-Free Systems
#
# Identify pairs with significant changes upon peptide binding (attractive/repulsive).
################################################################################

print("\n[4/8] Comparing peptide-loaded vs peptide-free systems...")
pairs = df_cons_pairs['resnum_chain12'].unique()
print(f"Analyzing {len(pairs)} pairs for peptide-binding effects...")

pairs_only_in_pmhc = list()
sig_aff_pairs_loaded = list()
p_values_loaded = list()

# Define the theoretical frame columns
frame_cols = list(map(str,np.arange(0, 251, 1)))
available_frame_cols = [col for col in frame_cols if col in df_intEn.columns]

for pair in tqdm.tqdm(pairs, desc="Testing peptide effects"):
    df_intEn_pair = df_cons_pairs[df_cons_pairs['resnum_chain12']==pair].copy()

    unique_systems = df_intEn_pair['system'].unique()
    unique_alleles = df_intEn_pair['allele'].unique()

    # Group systems by allele
    for allele in unique_alleles:
        df_intEn_pair_allele = df_intEn_pair[df_intEn_pair['allele']==allele].copy()
        unique_systems_allele = df_intEn_pair_allele['system'].unique()
        
        # Identify peptide-free system for this allele (system == allele, no peptide suffix)
        peptide_free_system = allele
        peptide_loaded_systems = [sys for sys in unique_systems_allele if sys != allele and sys.startswith(allele)]
        
        if peptide_free_system not in unique_systems_allele:
            # No peptide-free system for this allele, register peptide-loaded systems
            for sys in peptide_loaded_systems:
                pairs_only_in_pmhc.append([allele, sys, pair])
            continue
        
        if len(peptide_loaded_systems) == 0:
            # No peptide-loaded systems for this allele
            continue
            
        # Compare each peptide-loaded system with peptide-free system within the same allele
        for peptide_loaded_system in peptide_loaded_systems:
            df_pair_loaded = df_intEn_pair_allele[df_intEn_pair_allele['system']==peptide_loaded_system].copy()
            df_pair_free = df_intEn_pair_allele[df_intEn_pair_allele['system']==peptide_free_system].copy()
            
            if len(df_pair_loaded) == 0 or len(df_pair_free) == 0:
                continue

            val_arrays_loaded = list()
            val_arrays_free = list()

            for row in df_pair_loaded.iterrows():
                val_arrays_loaded.append(np.array(row[1][available_frame_cols]))

            val_arrays_loaded = np.hstack(val_arrays_loaded)

            for row in df_pair_free.iterrows():
                val_arrays_free.append(np.array(row[1][available_frame_cols]))

            val_arrays_free = np.hstack(val_arrays_free)

            # Calculate mean difference
            mean_diff = np.mean(val_arrays_loaded) - np.mean(val_arrays_free)
            if np.abs(mean_diff) < 2:
                # Magnitude-wise insignificant difference
                continue

            if mean_diff < 0:
                change_type = 'attractive'
            else:
                change_type = 'repulsive'

            # Perform statistical test
            try:
                t_statistic, p_value = scipy.stats.ttest_ind(val_arrays_loaded.astype('float'), 
                                                           val_arrays_free.astype('float'))
                if p_value < 0.05:
                    # Extract peptide from system name
                    peptide = peptide_loaded_system.replace(allele + '_', '') if '_' in peptide_loaded_system else peptide_loaded_system.replace(allele, '')
                    
                    sig_aff_pairs_loaded.append([allele, peptide, peptide_loaded_system, pair, change_type])
                    p_values_loaded.append(p_value)
                    
                    # Optional: Create visualization
                    # print(f'Significantly affected pair: {pair} in {allele} with peptide {peptide} (p-value: {p_value:.4f})')
                    # sns.kdeplot(val_arrays_loaded, color='blue', label=f'{allele}_{peptide}')
                    # sns.kdeplot(val_arrays_free, color='grey', label=f'{allele}_free')
                    # plt.legend()
                    # plt.title(f'Pair {pair} in {allele}: {change_type} change')
                    # plt.show()
                    
            except Exception as e:
                print(f"Error with pair {pair} in {allele} comparing {peptide_loaded_system} vs {peptide_free_system}: {e}")
                continue

df_sig_aff_pairs_loaded = pd.DataFrame(sig_aff_pairs_loaded, columns=['allele', 'peptide', 'system', 'pair', 'change_type'])
df_sig_aff_pairs_loaded['p_value'] = p_values_loaded

# Also create a DataFrame for pairs only in pMHC
df_pairs_only_in_pmhc = pd.DataFrame(pairs_only_in_pmhc, columns=['allele', 'system', 'pair'])

print(f"\nResults:")
print(f"  - {len(df_sig_aff_pairs_loaded)} significantly affected pairs")
print(f"  - {len(df_pairs_only_in_pmhc)} pairs only in peptide-loaded systems")
print(f"\nBreakdown by allele and change type:")
print(df_sig_aff_pairs_loaded.groupby(['allele', 'change_type']).size())
print(f"\nBreakdown by allele and peptide:")
print(df_sig_aff_pairs_loaded.groupby(['allele', 'peptide']).size())

print("\nSaving results...")
df_sig_aff_pairs_loaded.to_csv('df_sig_aff_pairs_loaded_2025_07_10.csv', index=False)
df_pairs_only_in_pmhc.to_csv('df_pairs_only_in_pmhc_2025_07_10.csv', index=False)
print("✓ Saved df_sig_aff_pairs_loaded_2025_07_10.csv")
print("✓ Saved df_pairs_only_in_pmhc_2025_07_10.csv")


################################################################################
# ## 1.5 Compare B*39:01/06 Subtypes with B*38:01
#
# Analyze differences between B*39 subtypes and B*38:01 for the same peptides.
################################################################################

print("\n[5/8] Comparing B*39 subtypes with B*38:01...")
pairs = df_cons_pairs['resnum_chain12'].unique()
print(f"Analyzing {len(pairs)} pairs for allele-specific differences...")

sig_aff_pairs_b39 = list()
p_values_b39 = list()

# Define the theoretical frame columns
frame_cols = list(map(str,np.arange(0, 251, 1)))
available_frame_cols = [col for col in frame_cols if col in df_intEn.columns]

for pair in tqdm.tqdm(pairs, desc="Testing allele differences"):
    df_intEn_pair = df_cons_pairs[df_cons_pairs['resnum_chain12']==pair].copy()

    unique_systems = df_intEn_pair['system'].unique()
    unique_alleles = df_intEn_pair['allele'].unique()
    unique_peptides = df_intEn_pair['peptide'].unique()
    # Only take strings in unique_peptides
    unique_peptides = [peptide for peptide in unique_peptides if isinstance(peptide, str) and peptide != 'nan']

    # Group systems by peptide
    for peptide in unique_peptides:
        df_intEn_pair_peptide = df_intEn_pair[df_intEn_pair['peptide']==peptide].copy()
        unique_systems_peptide = df_intEn_pair_peptide['system'].unique()
            
        # Compare each peptide-loaded system with peptide-free system within the same allele
        for peptide_loaded_system in unique_systems_peptide:
            df_pair_b39 = df_intEn_pair_peptide[df_intEn_pair_peptide['system']==peptide_loaded_system].copy()
            df_pair_b38 = df_intEn_pair_peptide[df_intEn_pair_peptide['allele']=='B3801'].copy()

            allele_b39 = df_pair_b39['allele'].unique()[0]

            if len(df_pair_b39) == 0 or len(df_pair_b38) == 0:
                continue

            val_arrays_b39 = list()
            val_arrays_b38 = list()

            for row in df_pair_b39.iterrows():
                val_arrays_b39.append(np.array(row[1][available_frame_cols]))

            val_arrays_b39 = np.hstack(val_arrays_b39)

            for row in df_pair_b38.iterrows():
                val_arrays_b38.append(np.array(row[1][available_frame_cols]))

            val_arrays_b38 = np.hstack(val_arrays_b38)

            # Calculate mean difference
            mean_diff = np.mean(val_arrays_b39) - np.mean(val_arrays_b38)
            if np.abs(mean_diff) < 2:
                # Magnitude-wise insignificant difference
                continue

            if mean_diff < 0:
                change_type = 'attractive'
            else:
                change_type = 'repulsive'

            # Perform statistical test
            try:
                t_statistic, p_value = scipy.stats.ttest_ind(val_arrays_b39.astype('float'), 
                                                           val_arrays_b38.astype('float'))
                if p_value < 0.05:
                    
                    sig_aff_pairs_b39.append([allele_b39, peptide, peptide_loaded_system, pair, change_type])
                    p_values_b39.append(p_value)
                    
                    # Optional: Create visualization
                    # print(f'Significantly affected pair: {pair} in {allele} with peptide {peptide} (p-value: {p_value:.4f})')
                    # sns.kdeplot(val_arrays_loaded, color='blue', label=f'{allele}_{peptide}')
                    # sns.kdeplot(val_arrays_free, color='grey', label=f'{allele}_free')
                    # plt.legend()
                    # plt.title(f'Pair {pair} in {allele}: {change_type} change')
                    # plt.show()
                    
            except Exception as e:
                print(f"Error with pair {pair} in {allele} comparing {peptide_loaded_system} vs B3801: {e}")
                continue

df_sig_aff_pairs_b39 = pd.DataFrame(sig_aff_pairs_b39, columns=['allele', 'peptide', 'system', 'pair', 'change_type'])
df_sig_aff_pairs_b39['p_value'] = p_values_b39

print(f"\nResults:")
print(f"  - {len(df_sig_aff_pairs_b39)} significantly affected pairs between B*39 and B*38:01")
print(f"\nBreakdown by allele and change type:")
print(df_sig_aff_pairs_b39.groupby(['allele', 'change_type']).size())
print(f"\nBreakdown by peptide and change type:")
print(df_sig_aff_pairs_b39.groupby(['peptide', 'change_type']).size())

print("\nSaving results...")
df_sig_aff_pairs_b39.to_csv('df_sig_aff_pairs_b39_2025_07_10.csv', index=False)
print("✓ Saved df_sig_aff_pairs_b39_2025_07_10.csv")


################################################################################
# ## 1.6 Construct Protein Energy Networks
#
# Implement Ribeiro-Ortiz method for building protein energy networks from interaction data.
################################################################################

print("\n[6/8] Constructing protein energy networks...")

def getRibeiroOrtizNetwork(pdb, df_intEn, includeCovalents=True,intEnCutoff=1,rmsdEn=5, startFrame=0):
	# Get the number of residues
	sys = parsePDB(pdb)
	sys_sel1 = sys.select("name CA")
	sys_sel2 = sys.select("name CA")
	numResidues_sel1 = len(np.unique(sys_sel1.getResindices()))
	numResidues_sel2 = len(np.unique(sys_sel2.getResindices()))

	# Get the system that includes both selections
	sys_sels = sys.select("name CA")
	resIndices_sels = np.unique(sys_sels.getResindices())
	numResidues = len(resIndices_sels)
	resIndices_sels_sel_str = ' '.join(list(map(str,resIndices_sels)))

	resNames = sys_sels.getResnames()
	resNums = sys_sels.getResnums()
	chains = sys_sels.getChids()

	# NOPE; THIS GET UNIQS APPROACH IS NOT WORKING. WE SHALL DO CA BASED LISTING!
	# WHY DID I NOT THINK ABOUT THIS BEFORE???
	def get_uniqs(original_list):
		
		result_list = list()
		seen = set()
		
		for item in original_list:
			if item not in seen:
				result_list.append(item)  # Add the first occurrence
				seen.add(item)
			elif item != result_list[-1]:
				result_list.append(item)  # Add subsequent occurrences of the last unique item

		return result_list

	rname_rnum_ch = zip(list(map(str,resNames)),list(map(str,resNums)),chains)
	rname_rnum_ch = ['_'.join(el) for el in rname_rnum_ch]
	#rname_rnum_ch = list(get_uniqs(rname_rnum_ch))

	# Get number of frames in df_intEn
	numFrames = len(df_intEn.columns[1:])
	# Start a list to store networks
	nx_list = list()

	for m in tqdm.tqdm(range(startFrame,numFrames), desc="Building networks", leave=False):

		frame_col = df_intEn.columns[1+m]
		# Create the network, and add nodes representing all residues in the protein
		network = nx.Graph()

		for j in range(0,numResidues):
			#print('adding: ',j+1,rname_rnum_ch[j])
			network.add_node(j+1, label=rname_rnum_ch[j])

		# Create an average residue interaction matrix for this specific frame.
		resIntEnMat = np.zeros((numResidues,numResidues))

		for row in df_intEn.iterrows():
			pair = row[1]['pair']
			resindex_1 = int(pair.split('-')[0])
			resindex_2 = int(pair.split('-')[1])
			resIntEnMat[resindex_1,resindex_2] = row[1][frame_col]
			resIntEnMat[resindex_2,resindex_1] = row[1][frame_col]

		# Determine RMSD of interaction energies (needed later on)
		#resIntEnArray = [i for i in np.reshape(resIntEnMat,(1,numResidues**2))[0]]
		#rmsdIntEn = np.sqrt(np.mean((resIntEnArray - np.mean(resIntEnArray)) ** 2))

		# Construct an matrix to make edge weights later on according to Ribeiro et al. (2014)
		#X = 0.5*(1-(resIntEnMat-np.mean(resIntEnArray))/(5*rmsdIntEn))

		# Construct an edge weights matrix.
		resIntEnMatNegFavor = np.zeros(np.shape(resIntEnMat))
		for i in range(0,np.shape(resIntEnMat)[0]):
			for j in range(0,np.shape(resIntEnMat)[0]):
				if resIntEnMat[i,j] < 0:
					resIntEnMatNegFavor[i,j] = np.abs(resIntEnMat[i,j])
				else:
					resIntEnMatNegFavor[i,j] = 0

		X = np.abs(resIntEnMatNegFavor)/np.max(np.abs(resIntEnMatNegFavor))

		for i in range(0,numResidues_sel1):
			for j in range(0,numResidues_sel2):

				if X[i,j] > 0.99:
					X[i,j] = 0.99

				if X[i,j] < 0:
					X[i,j] = 0
					print('alert: negative weight, reassigning to zero.')

		# If covalent bonds are also requested, then
		# Connect covalently bound residues with edge weight of 0.99 and an edge distance of 1/weight.
		if includeCovalents:
			for i in range(0,numResidues-1):
				resindex1 = resIndices_sels[i]
				resindex2 = resIndices_sels[i+1]
				res1 = sys.select('resindex %i' % resindex1)
				res2 = sys.select('resindex %i' % resindex2)
			
				# Only connect if the two residues are in the same chain
				# Weights as is
				if (res1.getChids()[0] == res2.getChids()[0]) and (res1.getSegindices()[0] == res2.getSegindices()[1]):
					network.add_edge(i+1,i+2,weight=X[i,i+1],distance=1-float(X[i,i+1]))


				# Weights with -log
				#if (res1.getChids()[0] == res2.getChids()[0]) and (res1.getSegindices()[0] == res2.getSegindices()[1]):
				#	network.add_edge(i+1,i+2,{'distance':-np.log(0.99)})

		for i in range(0,numResidues_sel1):
			for j in range(0,numResidues_sel2):

				if not includeCovalents:
				# Check again for covalent connection. If we are iterating over a residue pair covalently connected, then skip
					if abs(i-j) == 1:
						continue

			# Check whether edges exist between these residues.
				if not network.has_edge(i+1,j+1):

					# Check whether the mean interaction energy between the two residues is above the cutoff value.
					# If yes, continue.
					if abs(float(resIntEnMat[i,j])) >= float(abs(intEnCutoff)):

						# Add an edge between the two residues. Specify distance according to Ribeiro et al. (2014)
						# Also, consider edge weights lower than 0.01 disconnected
						# (again, Ribeiro et al. 2014)
						if X[i,j] < 0.01:
							continue

						# Connect the two residues with edge weight as calculated above and an edge distance of 1/weight.
						# weights as is
						network.add_edge(i+1,j+1,weight=X[i,j],distance=1-float(X[i,j]))
						# weights as -log
						#network.add_edge(i+1,j+1,{'distance':-np.log(X[i,j])})

		nx_list.append(network)

	return nx_list

def getBCs(dry_pdb_file,df_2network):
    nx_list = getRibeiroOrtizNetwork(dry_pdb_file,df_2network)
    df_bcs = [pd.DataFrame(nx.betweenness_centrality(nx_list[i]).items(),columns=['Residue','BC']) for i in range(0,len(nx_list))]
    i = 0
    for df in df_bcs:
        df['Frame'] = i
        nw = nx_list[i]
        node_labels = [nw.nodes[j]['label'] for j in df['Residue'].values]
        df['Label'] = node_labels
        i += 1
    df_bc = pd.concat(df_bcs,axis=0)
    return nx_list, df_bc

systems = df_intEn.system.unique()
numRepeats = 4
frame_cols = map(str,np.arange(0,1000,1))
available_frame_cols = [col for col in frame_cols if col in df_intEn.columns]

print(f"Processing {len(systems)} systems with up to {numRepeats} replicates each...")

df_bc_list = list()
intEn_cols = list(map(str,np.arange(0,880,10)))
for system in tqdm.tqdm(systems, desc="Computing betweenness centrality"):
    for i in range(1,numRepeats+1):
        df_2network = df_intEn.query(f"system == '{system}' & repeat == {i}")
        if len(df_2network) == 0:
            continue
        dry_pdb_file = df_2network['dry_pdb_file'].values[0]
        nx_list,df = getBCs(dry_pdb_file, df_2network[['pair']+available_frame_cols])
        df['system'] = system
        df['repeat'] = i
        df_bc_list.append(df)

print(f"Computed betweenness centrality for {len(df_bc_list)} system-replicate combinations")

df_bc = pd.concat(df_bc_list,axis=0,ignore_index=True)


################################################################################
# Compute betweenness centrality values for all residues across simulation frames.
################################################################################

df_bc['Residue_id'] = df_bc['Label'].str[4:]

print("Filtering equilibrated frames (Frame > 20)...")
df_bc_equil = df_bc.query('Frame > 20')
print(f"Retained {len(df_bc_equil)} BC measurements from equilibrated frames")

print("Saving to df_bc_equil_2025_07_10.csv...")
df_bc_equil.to_csv('df_bc_equil_2025_07_10.csv',index=False)
print("✓ Saved df_bc_equil_2025_07_10.csv")


################################################################################
# ## 1.7 Identify Residues with Consistent Betweenness Centrality
#
# Find residues maintaining stable BC values across simulation replicates (ANOVA p > 0.05).
################################################################################

print("\n[7/8] Identifying residues with consistent betweenness centrality...")
df_bc_equil.head()

systems = df_bc_equil.system.unique()
numRepeats = 4
print(f"Analyzing {len(systems)} systems for consistent residue BC...")

cons_resids_pdbs = list()
p_values = list()
for system in tqdm.tqdm(systems, desc="Analyzing systems"):
    df_bc_eq_pdb = df_bc_equil.query(f"system == '{system}'")

    resids = df_bc_eq_pdb.Residue_id.unique()

    for res in tqdm.tqdm(resids, desc=f"Testing residues in {system}", leave=False):
        df_res = df_bc_eq_pdb.query(f"Residue_id == '{res}'")

        # Get all available repeats for this residue
        available_repeats = df_res['repeat'].unique()

        # Need at least 2 repeats for ANOVA
        if len(available_repeats) < 2:
            continue

        # Collect BC values for each repeat
        bc_arrays = []
        for rep in available_repeats:
            bc_vals = df_res.query(f'repeat == {rep}')['BC'].values
            if len(bc_vals) > 0:  # Only include non-empty arrays
                bc_arrays.append(bc_vals)

        # Need at least 2 groups for ANOVA
        if len(bc_arrays) < 2:
            continue

        f_statistic, p_value = scipy.stats.f_oneway(*bc_arrays)
        if p_value > 0.05:
            cons_resids_pdbs.append([res,system])
            p_values.append(p_value)
            #print(f"Consistent residue: {lbl} in {pdb} simulations!")


################################################################################
# Find residues with consistent BC across replicates.
################################################################################

df_cons_resids_pdbs = pd.DataFrame(cons_resids_pdbs,columns=['Residue_id','system'])
df_cons_resids_pdbs['p_value'] = p_values

print(f"Found {len(cons_resids_pdbs)} residue-system pairs with consistent BC (p > 0.05)")

df_bc_equil_cons_resids = pd.merge(df_bc_equil,df_cons_resids_pdbs,how='outer',on=['Residue_id','system'])

print("Saving to df_bc_equil_cons_resids_2025_07_10.csv...")
df_bc_equil_cons_resids.to_csv('df_bc_equil_cons_resids_2025_07_10.csv',index=False)
print("✓ Saved df_bc_equil_cons_resids_2025_07_10.csv")


################################################################################
# ## 1.8 Identify Residues with Significant BC Changes
#
# Find residues whose network centrality changes significantly between peptide-binding conditions.
################################################################################

print("\n[8/8] Identifying residues with significant BC changes...")

################################################################################
# Assign peptide labels to systems for analysis.
################################################################################

print("Assigning peptide labels to systems...")
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3801_A','peptide'] = 'A'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3801_AL','peptide'] = 'AL'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3801_ALL','peptide'] = 'ALL'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3906_A','peptide'] = 'A'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3906_AL','peptide'] = 'AL'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3906_ALL','peptide'] = 'ALL'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3901_A','peptide'] = 'A'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3901_AL','peptide'] = 'AL'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3901_ALL','peptide'] = 'ALL'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3801','peptide'] = 'no peptide'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3906','peptide'] = 'no peptide'
df_bc_equil_cons_resids.loc[df_bc_equil_cons_resids["system"] == 'B3901','peptide'] = 'no peptide'

# Add allele column based on system name
df_bc_equil_cons_resids['allele'] = df_bc_equil_cons_resids['system'].apply(lambda x: x.split('_')[0] if '_' in x else x)

sig_aff_resids = list()
alleles = ['B3801', 'B3901', 'B3906']
peptides = ['A', 'AL', 'ALL']

unique_resids = df_cons_resids_pdbs['Residue_id'].unique()
print(f"Analyzing {len(unique_resids)} residues for significant BC changes...")

for resid in tqdm.tqdm(unique_resids, desc="Testing BC changes"):
    df_resid = df_bc_equil_cons_resids.query(f"Residue_id == '{resid}'")
    systems = df_resid.system.unique()

    if len(systems) == 1:
        continue
    
    # Build dict of BC values per system
    bc_resid_systems = dict()
    for system in systems:
        df_resid_system = df_resid.query(f"system == '{system}'")
        bc_resid_system = df_resid_system['BC'].values
        bc_resid_systems[system] = bc_resid_system

    # ANOVA across all systems for this residue
    f_statistic, p_value = scipy.stats.f_oneway(*list(bc_resid_systems.values()))
    
    if p_value < 0.05:
        result_row = {'Residue_id': resid, 'p_value': p_value}
        
        # For each allele, compare peptide-loaded to no-peptide reference
        for allele in alleles:
            ref_system = allele  # no peptide reference (e.g., 'B3801')
            
            if ref_system in bc_resid_systems:
                median_bc_ref = np.median(bc_resid_systems[ref_system])
                result_row[f'median_bc_{allele}_nopep'] = median_bc_ref
                
                for pep in peptides:
                    pep_system = f'{allele}_{pep}'  # e.g., 'B3801_A'
                    
                    if pep_system in bc_resid_systems:
                        median_bc_pep = np.median(bc_resid_systems[pep_system])
                        diff_bc = median_bc_pep - median_bc_ref
                        pt_diff_bc = diff_bc / np.abs(median_bc_pep) if median_bc_pep != 0 else None
                        
                        result_row[f'diff_bc_{allele}_{pep}'] = diff_bc
                        result_row[f'pt_diff_bc_{allele}_{pep}'] = pt_diff_bc
                    else:
                        result_row[f'diff_bc_{allele}_{pep}'] = None
                        result_row[f'pt_diff_bc_{allele}_{pep}'] = None
            else:
                result_row[f'median_bc_{allele}_nopep'] = None
                for pep in peptides:
                    result_row[f'diff_bc_{allele}_{pep}'] = None
                    result_row[f'pt_diff_bc_{allele}_{pep}'] = None
        
        sig_aff_resids.append(result_row)


################################################################################
# Compute differential BC between peptide-loaded and peptide-free conditions.
################################################################################

print(f"\nFound {len(sig_aff_resids)} residues with significant BC changes")

# Initialize common variables and load all preprocessed datasets
# sig_aff_resids is now a list of dicts, so DataFrame can infer columns automatically
df_bc_sig_aff_resids = pd.DataFrame(sig_aff_resids)

print("Saving to df_bc_sigaff_resids_2026_01_27.csv...")
df_bc_sig_aff_resids.to_csv('df_bc_sigaff_resids_2026_01_27.csv', index=False)
print("✓ Saved df_bc_sigaff_resids_2026_01_27.csv")


print("\n" + "="*80)
print("PREPROCESSING COMPLETE")
print("="*80)
print("\nGenerated CSV files:")
print("  - intEnVdW_2025_07_10.csv")
print("  - df_cons_pairs_saved_2025_07_10.csv")
print("  - df_sig_aff_pairs_loaded_2025_07_10.csv")
print("  - df_pairs_only_in_pmhc_2025_07_10.csv")
print("  - df_sig_aff_pairs_b39_2025_07_10.csv")
print("  - df_bc_equil_2025_07_10.csv")
print("  - df_bc_equil_cons_resids_2025_07_10.csv")
print("  - df_bc_sigaff_resids_2026_01_27.csv")
print("\nNext step: Run 2_analysis.ipynb for visualization and analysis")