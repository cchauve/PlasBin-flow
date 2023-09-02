""" Functions to compute parameters defining seeds """

import os
import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import defaultdict

from gfa_fasta_utils import (
    read_GFA_len
)

from ground_truth import (
    read_ground_truth_file
)

def _sample_ctg(sample, ctg):
    return f'{sample}_{ctg}'

def _read_sample_pls_score(sample, pls_score_file):
    """
    Reads plasmidness score file for a sample
    Returns dictionary:
	    Key: Contig id as sample_ctg
	    Value: Contig plasmidness (float)
    """
    pls_score_df = pd.read_csv(pls_score_file, sep='\t', header=None)
    pls_score_dict = dict(zip(pls_score_df[0], pls_score_df[1]))
    return {
        _sample_ctg(sample, ctg): pls_score
        for ctg, pls_score in pls_score_dict.items()
    }

def get_pls_ctgs_list(sample, ground_truth_df):	
    """
    Reads ground truth file for a sample
    Returns dictionary:
	    Key: Plasmid
	    Value: List of contig ids of contigs belonging to plasmid
	    Each contig id of the form sample_ctg
    """
    sample_pls_ctgs = ground_truth_df.groupby('plasmid')['contig'].apply(list).to_dict()
    for pls in sample_pls_ctgs:
        sample_pls_ctgs[pls] = [
            _sample_ctg(sample, ctg)
            for ctg in sample_pls_ctgs[pls]
        ]
    return sample_pls_ctgs

def change_type_to_pls(sample, ground_truth_df):
    """
    Function to change mol_type to 'plasmid' for 
    contigs mapped to plasmids in the ground truth.

    Takes ground truth dataframe as input
    Returns dictionary:
	    Key: Contig id as sample_ctg
	    Value: 'plasmid'
    """
    return {
        _sample_ctg(sample, ctg): 'plasmid'
        for ctg in ground_truth_df['contig'].unique()
    }

#Classifying contigs as seeds according to gene density and length thresholds
def seed_by_param(PLS_SC_THR, LEN_THR, CTG_LEN, CTG_PLS_SC):
    """
    Function to classify contigs as seeds
    according to given threshold params (gene density and length)
    """
    return 1 if CTG_PLS_SC >= PLS_SC_THR and CTG_LEN >= LEN_THR else 0	

def pls_seeds_by_thresholds(PLS_SC_THRESHOLDS, LEN_THRESHOLDS, CTG_DETAILS, PLS_CTGS):
    """
    Computing number of seeds for every plasmid
    Returns dictionary:
	    Key: Plasmid ID
	    Value: Nested dictionary with
		    Keys: threshold combinations (len_thr_sc_thr) 
		    Values: number of contigs classified as seeds using the thresholds
    """
    seeds_dict = {}
    for pls in PLS_CTGS:
        seeds_dict[pls] = {}
        for sc_thr in PLS_SC_THRESHOLDS:
            for len_thr in LEN_THRESHOLDS:
                seeds_dict[pls][f'{len_thr}_{sc_thr}'] = 0
                for ctg_id in PLS_CTGS[pls]:
                    seed_eligibility = seed_by_param(
                        sc_thr,
                        len_thr,
                        CTG_DETAILS['length'][ctg_id],
                        CTG_DETAILS['pls_score'][ctg_id]
                    )
                    seeds_dict[pls][f'{len_thr}_{sc_thr}'] += seed_eligibility
    return seeds_dict	

def count_false_seeds(PLS_SC_THRESHOLDS, LEN_THRESHOLDS, CTG_DETAILS_DF):
    """
    Computing contigs incorrectly classified as seeds
    Returns dictionary
	    Key: Plasmidness score threshold,
	    Value: Nested dictionary 
		    Keys: length thresholds
		    Value: number of false seeds
    """
    false_seeds_dict = {}
    for sc_thr in PLS_SC_THRESHOLDS:
        false_seeds_dict[sc_thr] = {}
        for len_thr in LEN_THRESHOLDS:
            false_seeds_dict[sc_thr][len_thr] = len(
                CTG_DETAILS_DF[
                    (CTG_DETAILS_DF['mol_type'] == 'chromosome') & \
                    (CTG_DETAILS_DF['pls_score'] >= sc_thr) & \
                    (CTG_DETAILS_DF['length'] >= len_thr)
                ]
            )
    return false_seeds_dict	

def count_pls_with_seeds(PLS_SC_THRESHOLDS, LEN_THRESHOLDS, seeds_df):
    """
    Computing number of plasmids with and without seed contigs
    Returns dictionary:
	    Key: Plasmidness score threshold,
	    Value: Nested dictionary 
		    Keys: Length thresholds
		    Value: number of plasmids with at least one seed
    """
    pls_with_seeds_dict = {}
    for sc_thr in PLS_SC_THRESHOLDS:
        pls_with_seeds_dict[sc_thr] = {}
        for len_thr in LEN_THRESHOLDS:
            pls_with_seeds_dict[sc_thr][len_thr] = len(
                seeds_df[seeds_df[f'{len_thr}_{sc_thr}'] >= 1]
            )
    return pls_with_seeds_dict

def output_best_params(pls_with_seeds_df, false_seeds_df, out_file):
    """
    Choosing seed parameters: plasmidness score threshold (sc_thr) and contig length (len_thr) threshold)
    We wish to choose seed parameters such that
    Number of seeded plasmids (SP) is maximized and
    Number of false seeds (NPS) is minimized.
    So, our objective is SP-NPS.
    """
    obj_df = pls_with_seeds_df.subtract(false_seeds_df)
    max_obj = obj_df.to_numpy().max()       #Computing the maximum objective value
    #Listing all combinations with max objective value
    with open(out_file, "w") as out:
        for row_name, row in obj_df.iterrows():
            for col_name, val in row.items():
                if val >= max_obj:
                    out.write(f'{row_name}\t{col_name}\n')	

def compute_seeds_parameters_file(input_csv_file, out_file):
    """
    Reads input csv file and creates a dataframe of file addresses 
    """
    PATHS_DF = pd.read_csv(
        input_csv_file,
        sep=',',
        names=[
            'sample','assembly','pls_score','ground_truth'
        ]
    )
    
    # Dictionaries for storing contig details and list of contigs belonging to plasmids
    CTG_DETAILS = {'length': {}, 'pls_score': {}, 'mol_type': {}}
    PLS_CTGS = {}
    for index, row in PATHS_DF.iterrows():
        # Reading sample GFA file and storing contig lengths
        sample_len_dict = read_GFA_len(
            row['assembly'], gzipped=True, id_fun=lambda x: _sample_ctg(row['sample'], x)
        )
        CTG_DETAILS['length'].update(sample_len_dict)
        # Reading and storing sample plasmidness scores
        sample_ctgs_list = list(sample_len_dict.keys())
        CTG_DETAILS['pls_score'].update(dict.fromkeys(sample_ctgs_list, 0))
        CTG_DETAILS['pls_score'].update(_read_sample_pls_score(row['sample'], row['pls_score']))
        # Reading and storing sample ground truth
        ground_truth_df = read_ground_truth_file(row['ground_truth'])
        PLS_CTGS.update(get_pls_ctgs_list(row['sample'], ground_truth_df))
        # Setting molecule type (chromosome/plasmid) for contigs
        CTG_DETAILS['mol_type'].update(dict.fromkeys(sample_ctgs_list, 'chromosome'))
        CTG_DETAILS['mol_type'].update(change_type_to_pls(row['sample'], ground_truth_df))
    # Ranges of thresholds
    LEN_THRESHOLDS = np.arange(50,5001,50)
    PLS_SC_THRESHOLDS = np.arange(1, 101, 1)/100
    seeds_dict = pls_seeds_by_thresholds(PLS_SC_THRESHOLDS, LEN_THRESHOLDS, CTG_DETAILS, PLS_CTGS)
    seeds_df = pd.DataFrame.from_dict(seeds_dict).T
    CTG_DETAILS_DF = pd.DataFrame.from_dict(CTG_DETAILS)
    false_seeds_dict = count_false_seeds(PLS_SC_THRESHOLDS, LEN_THRESHOLDS, CTG_DETAILS_DF)
    false_seeds_df = pd.DataFrame.from_dict(false_seeds_dict)
    pls_with_seeds_dict = count_pls_with_seeds(PLS_SC_THRESHOLDS, LEN_THRESHOLDS, seeds_df)
    pls_with_seeds_df = pd.DataFrame.from_dict(pls_with_seeds_dict)
    # Computing threshold pairs that maximize seeded plasmids and minimize false seeds and 
    # Writing to output file
    output_best_params(pls_with_seeds_df, false_seeds_df, out_file)
	






        
    
