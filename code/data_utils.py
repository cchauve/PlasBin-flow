""" Reading input data to process a sample with PlasBin-flow """

from __future__ import division
import math
from collections import defaultdict
import logging

from log_errors_utils import (
    process_exception,
    CustomException,
    check_lists_equality,
    check_lists_inclusion,
    check_num_fields,
    check_number_range
)
from gfa_fasta_utils import (
    read_GFA_len,
    read_GFA_normalized_coverage,
    read_GFA_links,
    GFA_TO_KEY,
    GFA_TO_ORIENT_KEY,
    GFA_FROM_ORIENT_KEY
)
from gc_content import (
    DEFAULT_GC_INTERVALS,
    read_gc_intervals_file,
    intervals_boundaries_to_str,
    read_gc_probabilities_file
)

"""
Contigs data structure:
Dictionary
Key: contig name
Values: Dictionary:
LEN_KEY: Length of the contig (int)
COV_KEY: Normalized coverage of the contig (float)
SCORE_KEY: Plasmid score (float)
SEED_KEY: Seed key (bool)
"""

LEN_KEY = 'length'
SCORE_KEY = 'plasmid_score'
COV_KEY = 'normalized_coverage'
SEED_KEY = 'seed'

UNICYCLER_TAG = 'unicycler'
SKESA_TAG = 'skesa'

# TAG of GFA file for recording normalized coverage
ASSEMBLER_COV_TAG = {
    UNICYCLER_TAG: 'dp',
    SKESA_TAG: None
}

# Default name for source and sink
DEFAULT_SOURCE = 'S'
DEFAULT_SINK = 'T'
# Head and tail string to define contig extremities
DEFAULT_HEAD_STR = 'h'
DEFAULT_TAIL_STR = 't'
# Default thresholds defining seeds used in the paper experiments
DEFAULT_SEED_LEN_THRESHOLD = 2650
DEFAULT_SEED_SCORE_THRESHOLD = 0.58
# Default plasmid score for short contigs or contigs with no score
DEFAULT_PLS_SCORE = 0.5
# Default minimum length under which contigs receive the default score
DEFAULT_MIN_CTG_LEN = 0

def read_pls_score_file(in_pls_score_file):
    """
    Reads a plasmid score file
    Args: 
        in_pls_score_file (str): path to a plasmid score file
    Returns: 
        Dictionary contig (str): plasmid score (float)
    Assumption:
        File existence and non-emptyness has been checked
    """
    pls_scores_dict = {}
    with open(in_pls_score_file) as in_file:
        for score_line in in_file.readlines():
            line = score_line.rstrip()
            line_split = line.split('\t')
            check_num_fields(line_split, 2)
            ctg_id = line_split[0]
            score = float(line_split[1])
            check_number_range(score, (0.0,1.0))
            pls_scores_dict[ctg_id] = score
    return pls_scores_dict

def read_ctgs_data(
    in_gfa_file, in_pls_score_file,
    assembler=UNICYCLER_TAG, gfa_gzipped=True,
    default_pls_score=DEFAULT_PLS_SCORE,
    min_ctg_len=DEFAULT_MIN_CTG_LEN    
):
    """
    Reads all data related to contigs
    Args:
        - in_gfa_file (str): path to a gzipped GFA file
        - in_pls_score_file (str): path to a plasmid score file
        - gfa_gzipped (bool): True if GFA file gzipped
        - assembler (str): tag of the used assembler
        - default_pls_score (float): default plasmid score
        - min_ctg_len (int): minimm length under which a contig 
                             receives the default plasmid score
    Returns:
        Dictionary described above without the SEEDS field
    """
    def _update_ctgs_dictionary(ctgs_data_dict, add_field_dict, add_key):
        for ctg_id,ctg_value in add_field_dict.items():
            ctgs_data_dict[ctg_id][add_key] = ctg_value

    ctgs_data_dict = defaultdict(dict)
    ctgs_len_dict = read_GFA_len(in_gfa_file, gzipped=gfa_gzipped)
    ctgs_cov_dict = read_GFA_normalized_coverage(
        in_gfa_file, gzipped=gfa_gzipped,
        cov_key=ASSEMBLER_COV_TAG[assembler]
    )
    _update_ctgs_dictionary(ctgs_data_dict, ctgs_len_dict, LEN_KEY)
    _update_ctgs_dictionary(ctgs_data_dict, ctgs_cov_dict, COV_KEY)

    ctgs_score_dict = read_pls_score_file(in_pls_score_file)
    check_lists_inclusion(
        ctgs_score_dict.keys(),
        ctgs_len_dict.keys(),
        msg=f'{in_pls_score_file} has contigs not present in {in_gfa_file}'
    )
    _update_ctgs_dictionary(ctgs_data_dict, ctgs_score_dict, SCORE_KEY)
    
    # If a contig is not in the plasmid scores file or is too short
    # it receives the default plasmid score
    for ctg_id,ctg_data in ctgs_data_dict.items():
        if SCORE_KEY not in ctg_data.keys() or ctg_data[LEN_KEY] < min_ctg_len:
            ctg_data[SCORE_KEY] = default_pls_score

    return ctgs_data_dict
    
def get_seeds(ctgs_data_dict, seed_len, seed_score,):
    """
    Computes the set of seed contigs and updates the contigs data
    Args:
        - ctgs_data_dict (Dictionary): see above
        - seed_len (int): length threshold defining seeds
        - seed_score float): plasmid score threshold defining seeds
    Returns:
        Set(contig names)
        Updates the fild SEED_KEY of ctgs_data_dict
    """
    seeds = set()
    for ctg_id,ctg_data in ctgs_data_dict.items():
        len_test = ctg_data[LEN_KEY] >= seed_len            
        score_test = ctg_data[SCORE_KEY] >= seed_score
        ctg_data[SEED_KEY] = len_test and score_test
        if ctg_data[SEED_KEY]:
            seeds.add(ctg_id)
    return seeds

def read_gc_data(gc_probabilities_file, gc_intervals_file, gfa_ctgs_list):
    """
    Reads the GC probabilities and intervals
    Computes the GC objective penalty for all contigs of a sample
    Args:
        gc_probabilities_file (str): path to GC probabilities file
        gc_intervals_file (str/None): path to GC intervals file
           if None, default intervals are used
        gfa_ctgs_list (List(str)) List of contig ids from GFA file
    Returns
        Dictionary contig id -> 
          Dictionary: interval string (str): GC probability for this interval
        Dictionary contig id -> 
          Dictionary: interval string (str): objective penalty for this interval
          Penalty = GC probability - max(GC probabilities for the contig)
    """
    gc_intervals_boundaries = (
        DEFAULT_GC_INTERVALS if gc_intervals_file is None 
        else read_gc_intervals_file(gc_intervals_file)
    )
    intervals_str_list = intervals_boundaries_to_str(gc_intervals_boundaries)
    num_intervals = len(intervals_str_list)

    __gc_probs_dict = read_gc_probabilities_file(
        gc_probabilities_file, num_intervals
    )
    gc_probs_dict,gc_pens_dict = {},{}
    for ctg_id,ctg_gc_probs in __gc_probs_dict.items():
        gc_probs_dict[ctg_id] = {
            intervals_str_list[i]: ctg_gc_probs[i]
            for i in range(len(intervals_str_list))
        }
        max_ctg_gc_prob = max(ctg_gc_probs)
        gc_pens_dict[ctg_id] = {
            int_str: ctg_gc_prob - max_ctg_gc_prob
            for int_str,ctg_gc_prob in gc_probs_dict[ctg_id].items()
        }
        
    check_lists_equality(
        gfa_ctgs_list,
        gc_probs_dict.keys(),
        msg=f'GFA file and {gc_probabilities_file} have inconsistent contig sets'
    )

    return gc_probs_dict,gc_pens_dict

def read_links_data(in_gfa_file, gfa_gzipped=True):
    """
    Reads a GFA file to extract edges and compute their capacity
    Args:
        - in_gfa_file (str): path to a GFA file
        - gfa_gzipped (bool): True if GFA file gzipped
    Returns:
        - List(
            (ctg1 (str),{DEFAULT_HEAD_STR,DEFAULT_TAIL_STR}),
            (ctg2 (str),{DEFAULT_HEAD_STR,DEFAULT_TAIL_STR})): list of edges
    """
    ctg_from_ext = {'+': DEFAULT_HEAD_STR, '-': DEFAULT_TAIL_STR}
    ctg_to_ext = {'+': DEFAULT_TAIL_STR, '-': DEFAULT_HEAD_STR}    
    links_list = []
    __links_dict = read_GFA_links(
        in_gfa_file,
        gzipped=gfa_gzipped
    )
    for ctg_id,ctg_links_from_list in __links_dict.items():
        for link in ctg_links_from_list:
            ext_from = (
                ctg_id,
                ctg_from_ext[link[GFA_FROM_ORIENT_KEY]]
            )
            ext_to = (
                link[GFA_TO_KEY],
                ctg_to_ext[link[GFA_TO_ORIENT_KEY]]
            )
            links_list.append((ext_from,ext_to))
    return links_list

def get_capacities(links_list, ctgs_data_dict):
    """
    Computes the capacity of all network edges
    Args:
        - links_list: see above
        - ctgs_data_dict: see above
    Returns:
        - Dictionary edge: capacity (float)
          edge is augmented by edges Source -> extremity and extremity -> Sink
    """
    capacities_dict = {}
    for link in links_list:
        ext_from,ext_to = link[0],link[1]
        ctg_from_id,ctg_to_id = ext_from[0],ext_to[0]
        capacity = min(
            ctgs_data_dict[ctg_from_id][COV_KEY],
            ctgs_data_dict[ctg_to_id][COV_KEY]
        )
        capacities_dict[(ext_from,ext_to)] = capacity
        capacities_dict[(ext_to,ext_from)] = capacity
    for ctg_id,ctg_data in ctgs_data_dict.items():
        ext_h,ext_t = (ctg_id, DEFAULT_HEAD_STR),(ctg_id, DEFAULT_TAIL_STR)
        capacity = ctg_data[COV_KEY]
        capacities_dict[(DEFAULT_SOURCE,ext_h)] = capacity
        capacities_dict[(DEFAULT_SOURCE,ext_t)] = capacity
        capacities_dict[(ext_h,DEFAULT_SINK)] = capacity
        capacities_dict[(ext_t,DEFAULT_SINK)] = capacity
    return capacities_dict

def log_data(ctgs_data_dict, links_list, in_gfa_file, in_pls_score_file):
    ctgs_list = ctgs_data_dict.keys()
    num_ctgs = len(ctgs_list)
    num_links = sum([len(links) for links in links_list])
    logging.info(
        f'DATA\tFile {in_gfa_file} contains {num_ctgs} contigs and {num_links} edges'
    )
    cov_list = [ctgs_data_dict[ctg_id][COV_KEY] for ctg_id in ctgs_list]
    logging.info(
        f'DATA\tFile {in_gfa_file} minimum coverage {min(cov_list)} maximum coverage {max(cov_list)}'
    )
    score_list = [ctgs_data_dict[ctg_id][SCORE_KEY] for ctg_id in ctgs_list]
    logging.info(f'DATA\tFile {in_pls_score_file} maximum plasmid score {max(score_list)}')
    num_seeds = len(
        [ctg_id for ctg_id in ctgs_list if ctgs_data_dict[ctg_id][SEED_KEY]]
    )
    if num_seeds == 0:
        logging.warning(f'DATA\tFile {in_gfa_file} has no seed')
    else:
        logging.info(f'DATA\tFile {in_gfa_file} has {num_seeds} seed(s)')
