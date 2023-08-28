''' Functions to manipulate seed contigs '''

import pandas as pd
import numpy as np

from gene_density import (
    compute_gene_density
)
from ground_truth import (
    read_ground_truth_file,
    GT_PLS_KEY,
    GT_CTG_KEY
)
from log_errors_utils import (
    process_exception
)

LENGTH_KEY = 'length'
GD_KEY = 'gd'
MOL_TYPE_KEY = 'mol_type'
MOL_TYPE_CHR = 'chromosome'
MOL_TYPE_PLS = 'plasmid'

SAMPLE_KEY = 'sample'

''' Reading input files '''

def _process_record(record):
    '''
    Computes sequence statistics

    Args:
        record (BioSeqIO.record): FASDTA record

    Returns:
        (Dictionary): LENGTH_KEY: length, GD_KEY: 0.0 (gene density), MOL_TYPE_KEY: MOL_TYPE_CHR
    '''
    return {
        LENGTH_KEY: len(str(record.seq)),
        GD_KEY: 0.0,
        MOL_TYPE_KEY: MOL_TYPE_CHR
    }

def _join_sample_ctg(sample, ctg):
    return f'{sample}_{ctg}'

''' Computing gene density for all contigs ''' 

def _compute_gene_density(sample, mappings_file, gfa_file, gfa_gzipped):
    '''
    Reading gene to contig mapping information
    Computes gene density for all contigs
    '''
    def _get_union(intervals):
        '''
        Takes the gene covering intervals for a contig and finds their union
        The length of the union is used to compute gene coverage
        '''
        intervals_union = []
        for _,begin,end in sorted(intervals):
            if intervals_union and intervals_union[-1][1] >= begin-1:
                intervals_union[-1][1] = max(intervals_union[-1][1],end)
            else:
                intervals_union.append([begin,end])
        return intervals_union

    def _compute_gd(intervals_union, ctg_len):
        '''
        Computes gene density using list of coverage intervals and contig length
        '''
        covered = 0
        for interval in intervals_union:
            covered += interval[1] - interval[0] + 1
        return covered / ctg_len
    
    mappings_df = read_blast_outfmt6_file(mappings_file)
    ctg_intervals = compute_blast_s_intervals(mappings_df)
    ctg_gd_dict = {}
    for ctg_id,intervals in ctg_intervals.items():
        intervals_union = _get_union(intervals)
        ctg_id_sample = _join_sample_ctg(sample, ctg_id)
        ctg_gd_dict[ctg_id_sample] = _compute_gd(intervals_union, ctg_len[ctg_id_sample])
    return ctg_gd_dict

''' Classifying contigs as seeds according to gene density and length thresholds '''
''' START TO CLEAN '''

def seed_by_param(gdt, lt, ctg_details):
    '''
    Function to classify contigs as seeds
    according to given threshold params (gene density and length)
    '''
    ctg_len, ctg_gd = ctg_details['length'], ctg_details['gd']
    return 1 if ctg_gd >= gdt and ctg_len >= lt else 0

def pls_seeds_by_thresholds(GD_THRESHOLDS, LEN_THRESHOLDS, all_ctgs_dict, all_pls_dict):
    '''
    Computing number of seeds for every plasmid
    Key: Plasmid ID
    Value: Nested dictionary with
    keys as theshold combinations (lt_gdt) and values as number of contigs classified as seeds using the thresholds
    '''
    seeds_dict = {}
    for pls in all_pls_dict:
        seeds_dict[pls] = {SAMPLE_KEY: all_pls_dict[pls][SAMPLE_KEY]}
        for gdt in GD_THRESHOLDS:
            for lt in LEN_THRESHOLDS:
                seeds_dict[pls][f'{lt}_{gdt}'] = 0
                for ctg_id in all_pls_dict[pls][CTG_LIST_KEY]:
                    seed_eligibility = seed_by_param(gdt, lt, all_ctgs_dict[ctg_id])
                    seeds_dict[pls][f'{lt}_{gdt}'] += seed_eligibility
    return seeds_dict

def count_false_seeds(GD_THRESHOLDS, LEN_THRESHOLDS, all_ctgs_df):
    '''
    Computing contigs incorrectly classified as seeds
    Key: Gene density threshold,
    Value: Nested dictionary with length thresholds as keys and number of false seeds as value
    '''
    false_seeds_dict = {}
    for gdt in GD_THRESHOLDS:
        false_seeds_dict[gdt] = {}
        for lt in LEN_THRESHOLDS:
            false_seeds_dict[gdt][lt] = len(
                all_ctgs_df[
                    (all_ctgs_df['type'] == 'chromosome') & \
                    (all_ctgs_df['gd'] >= gdt) & \
                    (all_ctgs_df['length'] >= lt)
                ]
            )
    return false_seeds_dict

def count_pls_with_seeds(GD_THRESHOLDS, LEN_THRESHOLDS, seeds_df):
    '''
    Computing number of plasmids with and without seed contigs
    Key: Gene density threshold,
    Value: Nested dictionary with length thresholds as keys and number of plasmids with seeds as values
    '''
    pls_with_seeds_dict = {}
    for gdt in GD_THRESHOLDS:
        pls_with_seeds_dict[gdt] = {}
        for lt in LEN_THRESHOLDS:
            pls_with_seeds_dict[gdt][lt] = len(
                seeds_df[seeds_df[f'{lt}_{gdt}'] >= 1]
            )
    return pls_with_seeds_dict

def output_best_params(pls_with_seeds_df, false_seeds_df, out_file):
    '''
    Choosing seed parameters
    We wish to choose seed parameters (gdt, lt) such that
    Number of seeded plasmids (SP) is maximized and
    Number of false seeds (NPS) is minimized.
    So, our objective is SP-NPS.
    '''
    obj_df = pls_with_seeds_df.subtract(false_seeds_df)
    
    max_obj = obj_df.to_numpy().max()       #Computing the maximum objective value
    
    #Listing all combinations with max objective value
    #best_params = set()
    with open(out_file, "w") as out:
        for row_name, row in obj_df.iterrows():
            for col_name, val in row.iteritems():
                if val >= max_obj:
                    #best_params.add(row_name, col_name)
                    out.write(f'{row_name}\t{col_name}\n')

''' END TO CLEAN '''
                    
def _read_ground_truth_file(all_ctgs_dict, all_pls_dict, ground_truth_file, sample):
    '''
    Reading ground truth file
    Dictionary with key as the plamid id
    Value as a dictionary of following attributes
    sample (str), ctg_list (list)
    '''
    try:
        ground_truth_df = read_ground_truth_file(ground_truth_file)
        for _,row in ground_truth_df.iterrows():
            pls, ctg = row[GT_PLS_KEY], row[GT_CTG_KEY]
            ctg_id = _join_sample_ctg(sample, ctg)
            all_ctgs_dict[ctg_id][MOL_TYPE_KEY] = MOL_TYPE_PLS
            if pls in all_pls_dict.keys():
                all_pls_dict[pls][CTG_LIST_KEY].append(ctg_id)
            else:
                all_pls_dict[pls] = {CTG_LIST_KEY: [ctg_id], SAMPLE_KEY: sample}
    except Exception as e:
        process_exception(
            f'Reading ground truth file {ground_truth_file}: {e}'
        )
    else:
        return all_ctgs_dict,all_pls_dict

def read_reference_data(input_csv_file):
    '''
    Takes csv file with addresses to assembly file,
    gene to contig mapping file (blast output) and ground truth file.
    Returns two dictionaries:
    1. Key: contig id, Value: nested dictionary with length, gd and contig source (plasmid/chromosome)
    2. Key: plasmid ids, Value: nested dictionary with sample name and list of contigs in the plasmid
    '''
    # Reading and storing input data for reference samples
    try:
        input_data_df = pd.read_csv(
            input_csv_file,
            sep = ',',
            header = None,
        )
    except Exception as e:
        process_exception(f'Reading CSV input file {input_csv_file}: {e}')
    
    all_ctgs_dict,all_pls_dict = {},{}
    for _, row in input_data_df.iterrows():
        sample, assembly, mappings, ground_truth = row[0], row[1], row[2], row[3]        
        all_ctgs_dict.update(
            _read_FASTA_file(assembly, sample)
        )
        ctg_len = {
            ctg_id: all_ctgs_dict[ctg_id][LENGTH_KEY]
            for ctg_id in all_ctgs_dict.keys()
        }
        ctg_gd_dict = compute_gene_density(sample, mappings, ctg_len)
        for ctg_id,ctg_gd in ctg_gd_dict.items():
            all_ctgs_dict[ctg_id][GD_KEY] = ctg_gd        
        all_ctgs_dict, all_pls_dict = _read_ground_truth_file(
            all_ctgs_dict, all_pls_dict, ground_truth, sample
        )
    return all_ctgs_dict,all_pls_dict


def compute_seeds_parameters_file(input_csv_file, out_file):
    #Reads reference data files and returns two dictionaries: one with contig details and another with plasmid details
    all_ctgs_dict, all_pls_dict = read_reference_data(input_csv_file)

    #Ranges of thresholds
    LEN_THRESHOLDS = np.arange(50,5001,50)
    GD_THRESHOLDS = np.arange(1, 101, 1)/100

    seeds_dict = pls_seeds_by_thresholds(GD_THRESHOLDS, LEN_THRESHOLDS, all_ctgs_dict, all_pls_dict)
    seeds_df = pd.DataFrame.from_dict(seeds_dict).T
    all_ctgs_df = pd.DataFrame.from_dict(all_ctgs_dict).T
    
    false_seeds_dict = count_false_seeds(GD_THRESHOLDS, LEN_THRESHOLDS, all_ctgs_df)
    false_seeds_df = pd.DataFrame.from_dict(false_seeds_dict)
    
    pls_with_seeds_dict = count_pls_with_seeds(GD_THRESHOLDS, LEN_THRESHOLDS, seeds_df)
    pls_with_seeds_df = pd.DataFrame.from_dict(pls_with_seeds_dict)
    
    output_best_params(pls_with_seeds_df, false_seeds_df, out_file)
