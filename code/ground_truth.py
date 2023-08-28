''' Functions to manipulate ground truth files '''

import pandas as pd

from mappings_utils import (
    read_blast_outfmt6_file,
    compute_blast_qs_intervals
)

GT_PLS_KEY = 'plasmid'
GT_CTG_KEY = 'contig'
GT_COV_KEY = 'coverage'
GT_PLSLEN_KEY = 'plasmid_len'
GT_CTGLEN_KEY = 'contig_len'

GT_COL_NAMES = [
    GT_PLS_KEY, GT_CTG_KEY, GT_COV_KEY, GT_PLSLEN_KEY, GT_CTGLEN_KEY
]
GT_COL_TYPES = {
    GT_PLS_KEY: str, GT_CTG_KEY: str, GT_COV_KEY: float,
    GT_PLSLEN_KEY: int, GT_CTGLEN_KEY: int
}

def _compute_true_positives(pls_intervals, ctg_len, cov_threshold):
    def _num_covered_positions(intervals):
        ''' 
        Computes the number of positions covered by a list of intervals 
        '''
        intervals.sort(key = lambda x: x[0])
        num_pos_covered,last_pos_covered = 0,0
        for start,end in intervals:
            if end > last_pos_covered:
                num_pos_covered += end - max(last_pos_covered + 1, start) + 1
                last_pos_covered = end
        return num_pos_covered
            
    TP_dict = {}
    for pls_id,ctgs_data in pls_intervals.items():
        for ctg_id,ctg_intervals in ctgs_data.items():
            ctg_coverage = _num_covered_positions(
                ctg_intervals
            )/ctg_len[ctg_id]
            if ctg_coverage >= cov_threshold:
                TP_dict[(pls_id,ctg_id)] = ctg_coverage
    return TP_dict

def _write_ground_truth_file(true_positive_dict, ctg_len, pls_len, ground_truth_file):
    with open(ground_truth_file, 'w') as out_file:
        for (pls_id,ctg_id),coverage in true_positive_dict.items():                    
            coverage_str = '{:.2f}'.format(coverage)
            out_file.write(
                f'{pls_id}\t{ctg_id}\t{coverage_str}\t{pls_len[pls_id]}\t{ctg_len[ctg_id]}\n'
            )

def read_ground_truth_file(ground_truth_file):
    return pd.read_csv(
        ground_truth_file,
        sep='\t',
        names=GT_COL_NAMES,
        dtype=GT_COL_TYPES
    )

def compute_ground_truth_file(
        out_dir, tmp_dir, sample, pls_mappings_file,
        ctg_len, pls_len, pid_threshold, cov_threshold,
        ground_truth_file):
    '''
    Computes the ground truth file for a sample

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       sample (str): sample -> sample id
       pls_mappings_file (str): path to file of mappings of contigs to true plasmids in BLAST format 6
       ctg_len (Dictionary): contig id -> contig length
       pls_len (Dictionary): plasmid id -> plasmid length
       pid_threshold (float): percent identity threshold
       cov_threshold (float): coverage threshold
       ground_truth_file (str): path to the ground truth file to write

    Returns:
      None, creates the file _ground_truth_file(out_dir, sample)
    '''
    pls_mappings_df = read_blast_outfmt6_file(
        pls_mappings_file, order_coordinates=True, min_pident=pid_threshold*100.0
    )
    pls_intervals = compute_blast_qs_intervals(pls_mappings_df)
    true_positive_dict = _compute_true_positives(pls_intervals, ctg_len, cov_threshold)
    _write_ground_truth_file(true_positive_dict, ctg_len, pls_len, ground_truth_file)

    
                
