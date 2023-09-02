""" Functions to manipulate mappings files """

import os
import glob
import pandas as pd
import logging
from log_errors_utils import (
    run_cmd,
    log_file
)

BLAST6_COL_NAMES = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
    'gapopen', 'qstart', 'qend', 'sstart', 'send',
    'evalue', 'bitscore'
]
BLAST6_COL_TYPES = {
    'qseqid': str, 'sseqid': str, 'pident': float, 'length': int, 'mismatch': int,
    'gapopen': int, 'qstart': int, 'qend': int, 'sstart': int, 'send': int,
    'evaue': float, 'bitscore': float
}

def run_blast6(query_file, db_file, mappings_file):
    """
    Blast query_file against db_file into mapping_file
    Args:
        - query_file (str): path to a FASTA file
        - db_file (str): path to a FASTA file
        - mappings_file (str): pth to the result file
    Returns:
        Creates mappings_file
    Assumptions:
       makeblastdb and blastn are in the default path
    """
    db_prefix = f'{db_file}.db'
    logging.info(f'ACTION\tcompute blast database for {db_file}: {db_prefix}')
    cmd_makeblastdb = [
        'makeblastdb',
        '-in', db_file,
        '-dbtype', 'nucl',
        '-out', db_prefix
    ]
    _ = run_cmd(cmd_makeblastdb)
    logging.info(f'ACTION\tmap {query_file} to {db_prefix}')        
    cmd_megablast = [
        'blastn', '-task', 'megablast',
        '-query', query_file,
        '-db', db_prefix,
        '-out', mappings_file,
        '-outfmt', '6'
    ]
    _ = run_cmd(cmd_megablast)
    log_file(mappings_file)
    db_files = glob.glob(f'{db_prefix}.n*')
    for file_to_clean in db_files:
        try:
            os.remove(file_to_clean)
        except Exception as e:
            process_exception('Removing BLAST file {file_to_clean}: (e)')

def _df_modify_col(in_df, col_name, modif_fun):
    """
    Modify in_df[col_name] inplace by applying modif_fun
    """
    for idx,row in in_df.iterrows():
        prev_val = row[col_name]
        new_val = modif_fun(prev_val)
        in_df.at[idx, col_name] = new_val

            
def read_blast_outfmt6_file(mappings_file, order_coordinates=True):
    """
    Reads a BLAST format 6 file

    Args:
        - mappings_file (str): BLAST format 6 output file
        - order_doordinates (bool): if True, reorder coordinates increasingly for each hit

    Returns:
        (DataFrame) with columns as in BLAST format 6 
        https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
        and start and end positions in increasing order
    """
    
    def _order_coordinates(in_df, start_col, end_col):
        """
        Ensures that the colums start_col,end_col are increasing
        """
        for idx,row in in_df.iterrows():
            start,end = row[start_col],row[end_col]
            if start>end:
                start,end=end,start
            in_df.at[idx,start_col] = start
            in_df.at[idx,end_col] = end
        
    mappings_df = pd.read_csv(
        mappings_file,
        sep='\t',
        names=BLAST6_COL_NAMES,
        dtype=BLAST6_COL_TYPES
    )
    if order_coordinates:
        _order_coordinates(mappings_df, 'qstart', 'qend')
        _order_coordinates(mappings_df, 'sstart', 'send')
    _df_modify_col(mappings_df, 'pident', modif_fun=lambda x: x/100.0)
    return mappings_df

def filter_blast_outfmt6(
        mappings_df, q_len_dict=None, s_len_dict=None,
        min_len=0, min_pident=0, min_q_cov=0.0, min_s_cov=0.0
):
    for mapping_index,mapping in mappings_df.iterrows():
        test_pident = mapping['pident'] < min_pident
        test_len = mapping['length'] < min_len
        q_id = mapping['qseqid']
        q_len = mapping['qend'] - mapping['qstart'] + 1
        test_q_cov = (q_len_dict) is not None and (q_len_dict[q_id]/q_len < min_q_cov)
        s_id = mapping['sseqid']
        s_len = mapping['send'] - mapping['sstart'] + 1
        test_s_cov = (s_len_dict) is not None and (s_len_dict[s_id]/s_len < min_s_cov)
        if test_pident or test_len or test_q_cov or test_s_cov:
            mappings_df.drop(mapping_index, inplace=True)
        
def compute_blast_qs_intervals(mappings_df):
    """
    From a mappings dataframe computes the intervals of each subject covered by each query

    Args: 
        mappings_df (DataFrame): mappings dataframe

    Returns:
        (Dictionary): subject id -> (Dictionary) query id -> List((qstart,qend,sstart,ssend)) of intervals 
                      of the subject covered by the query
    """
    s_ids = sorted(mappings_df.sseqid.unique())
    q_ids = sorted(mappings_df.qseqid.unique())
    q_s_ids = [
        (q_id,s_id)
        for q_id in q_ids
        for s_id in s_ids
    ]
    qs_intervals = {
        s_id: {
            q_id: []
            for q_id in q_ids
        }
        for s_id in s_ids
    }
    for (q_id,s_id) in q_s_ids:
        q_s_hits = mappings_df.loc[
            mappings_df.qseqid == q_id
        ].loc[mappings_df.sseqid == s_id]
        for _,row in q_s_hits.iterrows():
            interval = (row['qstart'], row['qend'], row['sstart'], row['send'])
            qs_intervals[s_id][q_id].append(interval)
    return qs_intervals

def compute_blast_s_intervals(mappings_df):
    """
    From a mappings dataframe computes the intervals of each subject covered by all queries

    Args: 
        mappings_df (DataFrame): mappings dataframe

    Returns:
        (Dictionary): subject id -> List((query,sstart,ssend)) of intervals 
    """
    s_ids = sorted(mappings_df.sseqid.unique())
    s_intervals = {s_id: [] for s_id in s_ids}
    for s_id in s_ids:
        s_hits = mappings_df.loc[mappings_df.sseqid == s_id]
        for _,row in s_hits.iterrows():
            interval = (row['qseqid'], row['sstart'], row['send'])
            s_intervals[s_id].append(interval)
    return s_intervals

def compute_blast_q_intervals(mappings_df):
    """
    From a mappings dataframe computes the intervals of each query covered by all subjects

    Args: 
        mappings_df (DataFrame): mappings dataframe

    Returns:
        (Dictionary): subject id -> List((subject,qstart,qend)) of intervals 
    """
    q_ids = sorted(mappings_df.qseqid.unique())
    q_intervals = {q_id: [] for q_id in q_ids}
    for q_id in q_ids:
        q_hits = mappings_df.loc[mappings_df.qseqid == q_id]
        for _,row in q_hits.iterrows():
            interval = (row['sseqid'], row['qstart'], row['qend'])
            q_intervals[q_id].append(interval)
    return q_intervals
