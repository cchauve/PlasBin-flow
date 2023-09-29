""" Functions to compute parameters defining seeds """

import os
import pandas as pd
import numpy as np
from collections import defaultdict
import logging

from gfa_fasta_utils import (
    read_GFA_len
)

from ground_truth import (
    read_ground_truth_file,
    GT_PLS_KEY,
    GT_CTG_KEY
)

from data_utils import (
    read_pls_score_file
)

# Should be added as parameters to seeds computation
SEEDS_DEFAULT_MIN_CTG_LEN=500
SEEDS_DEFAULT_MAX_CTG_LEN=5000
SEEDS_DEFAULT_CTG_LEN_STEP=50
SEEDS_DEFAULT_MIN_PLS_SCORE=0.01
SEEDS_DEFAULT_MAX_PLS_SCORE=1.0
SEEDS_DEFAULT_PLS_SCORE_STEP=0.01

"""
Reading data
"""

"""
all_ctgs_dict:

The first main data structure used in this module is a nested
dictionary with three inner dictionaries, each indexed by pairs (sample,ctg) 
for recording respectively, for each contig
- length of the contig (key: LEN_KEY),
- plasmid score of the contig (key: PLS_SCORE_KEY),
- molecule type of the contig (key: MOL_TYPE_KEY)
"""
MOL_TYPE_KEY = 'mol_type'
MOL_TYPE_CHR = 'chromosome'
MOL_TYPE_PLS = 'plasmid'
LEN_KEY = 'length'
PLS_SCORE_KEY = 'pls_score'

"""
pls_ctgs_dict:

The second main data structure is a dictionary indexed by true plasmids and
listing the contigs from that plasmid in the form (sample,contig)
"""

def _create_pls_ctgs_list(sample, ground_truth_df):	
    """
    Creates the list of contigs for each true plasmid for sample
    Args:
        - sample (str): sample name
        - ground_truth_df (DataFrame): ground truth dataframe obtained 
          by ground_truth.read_ground_truth_file
    Returns dictionary:
	Key: plasmid
	Value: List of contigs belonging to plasmid in format (sample,ctg)
    """
    pls_ctgs_dict = ground_truth_df.groupby(
        GT_PLS_KEY
    )[GT_CTG_KEY].apply(list).to_dict()
    return {
        (sample,pls_id): [(sample,ctg_id) for ctg_id in pls_ctgs]
        for pls_id,pls_ctgs in pls_ctgs_dict.items()
    }

def _read_input_data(input_file):
    """
    Reads the seeds input file and returns dictionary on contigs and plasmids
    Args:
       - input_file: path to CSV file with fields 
         sample, gfa_file, pls_score_file, ground_truth_file
    Returns
       - Dictionary: 
         LEN_KEY:       Dictionary (sample,contig): length
         PLS_SCORE_KEY: Dictionary (sample,contig): plasmid score
         MOL_TYPE_KEY:  Dictionary (sample,contig): plasmid/chromsome
       - Dictionary (sample,plasmid): list (sample,contig) of contigs in plasmid ground truth
    """
    SAMPLE_KEY = 'sample'
    GFA_FILE_KEY = 'gfa_file'
    PS_FILE_KEY = 'pls_score_file'
    GT_FILE_KEY = 'ground_truth_file'
    input_df = pd.read_csv(
        input_file,
        sep=',',
        names=[SAMPLE_KEY,GFA_FILE_KEY,PS_FILE_KEY,GT_FILE_KEY]
    )
    all_ctgs_dict = {
        LEN_KEY: {},  PLS_SCORE_KEY: {}, MOL_TYPE_KEY: {}
    }
    pls_ctgs_dict = {}
    for _,sample_row in input_df.iterrows():
        sample = sample_row[SAMPLE_KEY]
        # Reading all contigs length from sample GFA file
        sample_ctgs_len_dict = read_GFA_len(
            sample_row[GFA_FILE_KEY], gzipped=True
        )
        for ctg_id,ctg_len in sample_ctgs_len_dict.items():
            all_ctgs_dict[LEN_KEY][(sample,ctg_id)] = ctg_len
            all_ctgs_dict[PLS_SCORE_KEY][(sample,ctg_id)] = 0.0
            all_ctgs_dict[MOL_TYPE_KEY][(sample,ctg_id)] = MOL_TYPE_CHR
        # Updating contigs plasmid score from sample plasmid scores file
        sample_ctgs_pls_score_dict = read_pls_score_file(sample_row[PS_FILE_KEY])
        for ctg_id,ctg_pls_score in sample_ctgs_pls_score_dict.items():
            all_ctgs_dict[PLS_SCORE_KEY][(sample,ctg_id)] = ctg_pls_score
        # Reading sample ground truth
        sample_ground_truth_df = read_ground_truth_file(sample_row[GT_FILE_KEY])
        # Recording the contigs on true plasmid
        for _,gt_row in sample_ground_truth_df.iterrows():
            ctg_id = gt_row[GT_CTG_KEY]
            all_ctgs_dict[MOL_TYPE_KEY][(sample,ctg_id)] = MOL_TYPE_PLS
        # Updating the dictionary (sample,plasmid) -> list of contigs on plasmid
        sample_pls_ctgs_list_dict = _create_pls_ctgs_list(
            sample, sample_ground_truth_df
        )
        for (sample,pls_id),pls_ctgs_list in sample_pls_ctgs_list_dict.items():
            pls_ctgs_dict[(sample,pls_id)] = pls_ctgs_list
    return all_ctgs_dict,pls_ctgs_dict
        
"""
Defining seeds
"""

def count_false_seeds(thresholds, all_ctgs_dict):
    """
    Computing contigs incorrectly classified as seeds, that are
    contigs that pass both thrsholds but are assigned to chromosomes
    Args:
        - thresholds (List((int,float)): list of (length,score) thresholds
        - all_ctgs_dict (Dictionary): see all_ctgs_dict description above    
    Returns Dictionary 
       Key: pair of thresholds (int,float)
       Value: list of false seeds (sample,ctg)
    """
    all_ctgs = all_ctgs_dict[LEN_KEY].keys()
    false_seeds_dict = {}
    for (len_thr,pls_score_thr) in thresholds:
        false_seeds_dict[(len_thr,pls_score_thr)] = []
        for ctg_id in all_ctgs:
            test_chr = all_ctgs_dict[MOL_TYPE_KEY][ctg_id] == MOL_TYPE_CHR
            test_len = all_ctgs_dict[LEN_KEY][ctg_id] >= len_thr
            test_score = all_ctgs_dict[PLS_SCORE_KEY][ctg_id] >= pls_score_thr
            if test_chr and test_len and test_score:
                false_seeds_dict[(len_thr,pls_score_thr)].append(ctg_id)
    return false_seeds_dict

def count_pls_with_seed(thresholds, all_ctgs_dict, pls_ctgs_dict):
    """
    Computing number of plasmids with at least one seed contig
    Args:
        - thresholds (List((int,float)): list of (length,score) thresholds
        - all_ctgs_dict (Dictionary): see all_ctgs_dict description above 
        - pls_ctgs_dict (Dictionary): see pls_ctgs_dict description above 
    Returns 
       - Dictionary:
         Key: pair of thresholds (int,float)
         Value: list of (sample,plasmid) with at least one seed
       - Dictionary:
         Key: pair of thresholds (int,float)
         Value: list of (sample,plasmid) with at no seed
    """
    pls_with_seed_dict = {}
    pls_without_seed_dict = {}
    for (len_thr,pls_score_thr) in thresholds:
        pls_with_seed_dict[(len_thr,pls_score_thr)] = []
        pls_without_seed_dict[(len_thr,pls_score_thr)] = []
        for sample_pls_id,pls_ctgs_list in pls_ctgs_dict.items():
            add_pls_with_seed = False
            for ctg_id in pls_ctgs_list:
                ctg_len = all_ctgs_dict[LEN_KEY][ctg_id]
                ctg_pls_score = all_ctgs_dict[PLS_SCORE_KEY][ctg_id]
                if  ctg_pls_score >= pls_score_thr and ctg_len >= len_thr:
                    add_pls_with_seed = True
            if add_pls_with_seed:
                pls_with_seed_dict[(len_thr,pls_score_thr)].append(sample_pls_id)
            else:
                pls_without_seed_dict[(len_thr,pls_score_thr)].append(sample_pls_id)
    return pls_with_seed_dict,pls_without_seed_dict

def select_best_thresholds(thresholds, all_ctgs_dict, pls_ctgs_dict):
    """
    Choosing seed parameters: plasmid score threshold  and contig length threshold
    We wish to choose threshold pairs such that
    Number of seeded plasmids (SP) is maximized and
    Number of false seeds (NPS) is minimized.
    Our objective is to maximize SP-NPS.
    Args:
        - thresholds (List((int,float)): list of (length,score) thresholds
        - all_ctgs_dict (Dictionary): see all_ctgs_dict description above 
        - pls_ctgs_dict (Dictionary): see pls_ctgs_dict description above 
    Returns
        - List((int,float,int,int,int)): list of optimal threshold pairs (len,score)
          augmented by the number of plasmid with seed, the number of false seeds
          and the number of plasmids without seed
        - List((int,float,int,int,int)): all results
    """
    pls_with_seed_dict,pls_without_seed_dict = count_pls_with_seed(
        thresholds, all_ctgs_dict, pls_ctgs_dict
    )
    false_seeds_dict = count_false_seeds(thresholds, all_ctgs_dict)
    best_objective,best_thresholds,all_thresholds = None,[],[]
    for thr_pair in thresholds:
        pls_with_seed = pls_with_seed_dict[thr_pair]
        pls_without_seed = pls_without_seed_dict[thr_pair]
        false_seeds = false_seeds_dict[thr_pair]
        objective = len(pls_with_seed) - len(false_seeds)
        result = [
            thr_pair[0], thr_pair[1],
            len(pls_with_seed), len(false_seeds), len(pls_without_seed)
        ]
        if best_objective is None or objective > best_objective:
            best_thresholds = [result]
            best_objective = objective
        elif objective == best_objective:
            best_thresholds.append(result)
        all_thresholds.append(result)
    return best_thresholds,all_thresholds

""" I/O """
            
def _write_seeds_thresholds_file(thresholds, out_file_name):
    with open(out_file_name, "w") as out_file:
        out_file.write(
            f'contig_length\tplasmid_score\tseed_score\tplasmids_with_seed\tfalse_seeds\tplasmids_without_seed\n'
        )
        for [len_thr,score_thr,pls_with_seed,false_seeds,pls_without_seed] in thresholds:
            rounded_score_thr = round(score_thr,2)
            seed_score = pls_with_seed - false_seeds
            out_file.write(
                f'{len_thr}\t{rounded_score_thr}\t{seed_score}\t{pls_with_seed}\t{false_seeds}\t{pls_without_seed}\n'
            )

def read_seeds_thresholds_file(in_file_name):
    thresholds = []
    with open(in_file_name) as in_file:
        for line in in_file.readlines():
            line_split = line.restrip().split('\t')
            if line_split[0] != 'contig_length':
                thresholds.append(
                    [
                        int(line_split[0]),
                        float(line_split[1]),
                        int(line_split[2]),
                        int(line_split[3]),
                        int(line_split[4]),
                        int(line_split[5])
                    ]
                )
    return thresholds

""" Main """

def compute_optimal_seeds_parameters(
        input_file,
        min_ctg_len, max_ctg_len, ctg_len_step,
        min_pls_score, max_pls_score, pls_score_step
):
    """
    Compute optimal seeds parameters
    """
    all_ctgs_dict,pls_ctgs_dict = _read_input_data(input_file)
    combined_thresholds = [
        (len_threshold,pls_score_threshold)
        for len_threshold in np.arange(
                min_ctg_len, max_ctg_len+ctg_len_step, ctg_len_step
        )
        for pls_score_threshold in np.arange(
                min_pls_score, max_pls_score+pls_score_step, pls_score_step
        )
    ]
    best_thresholds,all_thresholds = select_best_thresholds(
        combined_thresholds, all_ctgs_dict, pls_ctgs_dict
    )
    return best_thresholds,all_thresholds

def compute_seeds_parameters_file(
        input_file, out_file_name,
        min_ctg_len=SEEDS_DEFAULT_MIN_CTG_LEN,
        max_ctg_len=SEEDS_DEFAULT_MAX_CTG_LEN,
        ctg_len_step=SEEDS_DEFAULT_CTG_LEN_STEP,
        min_pls_score=SEEDS_DEFAULT_MIN_PLS_SCORE,
        max_pls_score=SEEDS_DEFAULT_MAX_PLS_SCORE,
        pls_score_step=SEEDS_DEFAULT_PLS_SCORE_STEP
):
    best_thresholds,all_thresholds = compute_optimal_seeds_parameters(
        input_file,
        min_ctg_len, max_ctg_len, ctg_len_step,
        min_pls_score, max_pls_score, pls_score_step
    )
    if len(best_thresholds) == 0:
        logging.warning('No optimal seeds paramaters were found')
    else:
        logging.info(
            f'{len(best_thresholds)} optimal seeds paramaters were found'
        )
    _write_seeds_thresholds_file(best_thresholds, out_file_name)
    (out_file_name_prefix,out_file_name_ext) = os.path.splitext(out_file_name)
    out_file_name_all = os.path.normpath(
        f'{out_file_name_prefix}_all{out_file_name_ext}'
    )
    _write_seeds_thresholds_file(all_thresholds, out_file_name_all)
