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

from get_data import (
    read_pls_score_file
)

# Default thresholds defining seeds used in the paper experiments
DEFAULT_SEED_LEN_THRESHOLD = 2650
DEFAULST_SEED_SCORE_THRESHOLD = 0.58

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

def _read_sample_ctgs_len(sample, gfa_file):
    """
    Reads contigs name and length from a GFA file
    Args:
        - sample (str): sample name
        - gfa_file (str): path to a gzipped GFA file
    Returns dictionary:
        Key (sample,ctg) 
        Value: contig length (int)
    """
    return {
        (sample,ctg_id): ctg_len
        for ctg_id,ctg_len in read_GFA_len(gfa_file, gzipped=True).items()
    }

def _read_sample_pls_score(sample, pls_score_file):
    """
    Reads plasmid score file for a sample
    Args:
        - sample (str): sample name
        - pls_score_file (str): path to a plasmid score file
    Returns Dictionary 
        Key: (sample,ctg)
        Value: contig plasmid score (float) 
    """
    return {        
        (sample,ctg_id): pls_score
        for ctg_id,pls_score in read_pls_score_file(
                pls_score_file
        ).items()
    }

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
	Key: Plasmid
	Value: List of contigs belonging to plasmid in format (sample,ctg)
    """
    pls_ctgs_dict = ground_truth_df.groupby(
        GT_PLS_KEY
    )[GT_CTG_KEY].apply(list).to_dict()
    return {
        pls_id: [(sample,ctg_id) for ctg_id in pls_ctgs]
        for pls_id,pls_ctgs in pls_ctgs_dict.items()
    }

def _read_input_data(input_file):
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
        LEN_KEY: {},  PLS_SCORE_KEY: {},
        MOL_TYPE_KEY: defaultdict(lambda: MOL_TYPE_CHR)
    }
    pls_ctgs_dict = {}
    for _,sample_row in input_df.iterrows():
        all_ctgs_dict[PLS_SCORE_KEY].update(
            _read_sample_pls_score(
                sample_row[SAMPLE_KEY], sample_row[PS_FILE_KEY]
            )
        )
        all_ctgs_dict[LEN_KEY].update(
            _read_sample_ctgs_len(
                sample_row[SAMPLE_KEY], sample_row[GFA_FILE_KEY]
            )
        )
        sample_ground_truth_df = read_ground_truth_file(
            sample_row[GT_FILE_KEY]
        )
        for _,gt_row in sample_ground_truth_df.iterrows():
            ctg_key = (sample_row[SAMPLE_KEY], gt_row[GT_CTG_KEY])
            all_ctgs_dict[MOL_TYPE_KEY][ctg_key] = MOL_TYPE_PLS
        pls_ctgs_dict.update(
            _create_pls_ctgs_list(
                sample_row[SAMPLE_KEY], sample_ground_truth_df
            )
        )
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
       Value: number of false seeds
    """
    all_ctgs = all_ctgs_dict[LEN_KEY].keys()
    false_seeds_dict = defaultdict(int)
    for (len_threshold,score_threshold) in thresholds:
        for ctg_id in all_ctgs:
            test_chr = all_ctgs_dict[MOL_TYPE_KEY][ctg_id] == MOL_TYPE_CHR
            test_len = all_ctgs_dict[LEN_KEY][ctg_id] >= len_threshold
            test_score = all_ctgs_dict[PLS_SCORE_KEY][ctg_id] >= score_threshold
            if test_chr and test_len and test_score:
                false_seeds_dict[(len_threshold,score_threshold)] += 1
    return false_seeds_dict

def count_pls_with_seeds(thresholds, all_ctgs_dict, pls_ctgs_dict):
    """
    Computing number of plasmids with at least one seed contig
    Args:
        - thresholds (List((int,float)): list of (length,score) thresholds
        - all_ctgs_dict (Dictionary): see all_ctgs_dict description above 
        - pls_ctgs_dict (Dictionary): see pls_ctgs_dict description above 
    Returns Dictionary:
       Key: pair of thresholds (int,float)
       Value: number of plasmids with at least one seed
    """
    pls_with_seeds_dict = defaultdict(int)
    for (len_thr,pls_score_thr) in thresholds:
        for pls_ctgs_list in pls_ctgs_dict.values():
            for ctg_id in pls_ctgs_list:
                ctg_len = all_ctgs_dict[LEN_KEY][ctg_id]
                ctg_pls_score = all_ctgs_dict[PLS_SCORE_KEY][ctg_id]
                if  ctg_pls_score >= pls_score_thr and ctg_len >= len_thr:
                    pls_with_seeds_dict[(len_thr,pls_score_thr)] += 1
                    break
    return pls_with_seeds_dict

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
        List((int,float,int,int)): list of optimal threshold pairs (len,score)
        augmented by the number of plasmid with seeds and the number of false seeds
    """
    pls_with_seeds_dict = count_pls_with_seeds(
        thresholds, all_ctgs_dict, pls_ctgs_dict
    )
    false_seeds_dict = count_false_seeds(thresholds, all_ctgs_dict)
    best_objective,best_thresholds = None,[]
    for thr_pair in thresholds:
        pls_with_seeds = pls_with_seeds_dict[thr_pair]
        false_seeds = false_seeds_dict[thr_pair]
        objective = pls_with_seeds - false_seeds
        if best_objective is None or objective > best_objective:
            best_thresholds = [
                [thr_pair[0], thr_pair[1], pls_with_seeds, false_seeds]
            ]
            best_objective = objective
        elif objective == best_objective:
            best_thresholds.append(
                [thr_pair[0], thr_pair[1], pls_with_seeds, false_seeds]
            )
    return best_thresholds

""" I/O """
            
def _write_seeds_thresholds_file(thresholds, out_file_name):
    with open(out_file_name, "w") as out_file:
        for [len_thr,score_thr,pls_with_seeds,false_seeds] in thresholds:
            out_file.write(
                f'{len_thr}\t{score_thr}\t{pls_with_seeds}\t{false_seeds}\n'
            )

def read_seeds_thresholds_file(in_file_name):
    thresholds = []
    with open(in_file_name) as in_file:
        for line in in_file.readlines():
            line_split = line.restrip().split('\t')
            thresholds.append(
                [
                    int(line_split[0]),
                    float(line_split[1]),
                    int(line_split[2]),
                    int(line_split[3])
                ]
            )
    return thresholds

""" Main """

def compute_optimal_seeds_parameters(
        input_file, min_len, max_len, len_step, pls_score_step
):
    """
    Compute optimal seeds parameters
    """
    all_ctgs_dict,pls_ctgs_dict = _read_input_data(input_file)
    combined_thresholds = [
        (len_threshold,pls_score_threshold)
        for len_threshold in np.arange(min_len,max_len+1,len_step)
        for pls_score_threshold in np.arange(1, 101, pls_score_step)/100
    ]
    best_thresholds = select_best_thresholds(
        combined_thresholds, all_ctgs_dict, pls_ctgs_dict
    )
    return best_thresholds

def compute_seeds_parameters_file(
        input_file, out_file_name,
        min_len=50, max_len=5000, len_step=50, pls_score_step=1
):
    best_thresholds = compute_optimal_seeds_parameters(
        input_file, min_len, max_len, len_step, pls_score_step
    )
    if len(best_thresholds) == 0:
        logging.warning('No optimal seeds paramaters were found')
    else:
        logging.info(
            f'{len(best_thresholds)} optimal seeds paramaters were found'
        )
    _write_seeds_thresholds_file(best_thresholds, out_file_name)
