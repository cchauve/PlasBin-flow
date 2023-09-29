""" Functions to compute GC content intervals and probabilities """

import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.special as sc
import scipy.integrate as integrate
import numpy as np
from scipy.stats import binom
import math

from gfa_fasta_utils import (
    read_FASTA_ctgs,
    read_GFA_ctgs,
    GFA_SEQ_KEY
)
from log_errors_utils import (
    CustomException,
    process_exception,
    check_num_fields_eq,
    check_number_range,
    check_number_eq,
    check_number_lt
)

# Keys of dictionary associated to a FASTA record
GC_COUNT_KEY = 'gc_count'
GC_RATIO_KEY = 'gc_ratio'
LENGTH_KEY = 'length'
MOL_TYPE_KEY = 'mol_type'
MOL_TYPE_CHR = 'chromosome'
MOL_TYPE_PLS = 'plasmid'

DEFAULT_GC_INTERVALS = [0, 0.4,  0.45, 0.5, 0.55, 0.6, 1]
DEFAULT_N_GCINTS = 6

def _process_sequence(sequence):
    """
    Computes sequence statistics

    Args:
        sequence (str): input sequence
    
    Returns:
        (Dictionary): GC_COUNT_KEY -> GC content, GC_RATIO_KEY: GC content ratio, LENGTH_KEY: length
    """
    gc_count = sequence.count('G') + sequence.count('C')
    length = len(sequence)
    return {
        GC_COUNT_KEY: gc_count,
        GC_RATIO_KEY: gc_count / length if length > 0 else 0.0,
        LENGTH_KEY: length
    }

def _process_FASTA_record(record):
    """
    Parse a FASTA record

    Args:
        record (Bio.SeqIO.record): FASTA record

    Returns:
        (Dictionary): GC count (GC_COUNT_KEY), length (LENGTH_KEY), GC ratio (GC_RATIO_KEY)
    """
    return _process_sequence(str(record.seq))
    

def _process_GFA_ctg(attributes_dict):
    """
    Parse a GFA record

    Args:
        (Dictionary): attribute key: attribute value
        where attribute key is GFA_SEQ_KEY (contig sequence)

    Returns:
        (Dictionary): GC count (GC_COUNT_KEY), GC ratio (GC_RATIO_KEY), length (LENGTH_KEY)
    """
    return _process_sequence(attributes_dict[GFA_SEQ_KEY])

def _read_FASTA_files(fasta_file_paths, mol_type):
    """
    Read all plasmid/chromosome FASTA files, recording GC features

    Args:
        - fasta_file_paths (List(str)): list of paths to FASTA files to read
        - mol_type (str): MOL_TYPE_CHR, MOL_TYPE_PLS

    Returns:
        (Dictionary): sequence id ->
            (Dictionary) feature key -> feature
            feature key in {GC_COUNT_KEY, GC_RATIO_KEY, LENGTH_KEY, MOL_TYPE_KEY}
    """
    seq_dict = {}
    for fasta_file in fasta_file_paths:
        seq_dict.update(
            read_FASTA_ctgs(
                fasta_file, _process_FASTA_record, gzipped=True
            )
        )
        for seq_id in seq_dict.keys():
            seq_dict[seq_id][MOL_TYPE_KEY] = mol_type
    return seq_dict

def compute_gc_intervals_files(
        chr_paths_file, pls_paths_file, out_csv_file, out_png_file, out_intervals_file, n_gcints
):
    """
    Reads fasta sequences of plasmids and chromosomes in the reference database
    Computes and outputs details on gc content of molecules including
    violinplot of gc content distribution

    Args:
        - chr_paths_file (str): path to file containing the list of paths of chromosome FASTA files
        - pls_paths_file (str): path to file containing the list of paths of plasmids FASTA files
        - out_csv_file (str): written CSV file containing GC content of chromosomes and plasmids
        - out_png_file (str): written PNG file containing violin plot of GC content of chromosomes and plasmids
        - out_intervals_file (str): written text file containing the endpoints of GC content intervals in ascending order
        - n_gcints (int): number of GC content intervals between 0 and 1

    Returns:
        None, creates files out_csv_file, out_png_file, out_intervals_file
        out_csv_file and out_png_file record GC content and ratio features
        out_intervals_file records n_gcints equally spaced intervals between the min and mac GC content ratio
    """

    # Reading and storing input data for reference samples
    # seq_dict: sequence id -> {
    #  GC_COUNT_KEY -> GC content, GC_RATIO_Key: GC content ratio, LENGTH_KEY: length
    # }
    chr_paths = open(chr_paths_file,'r').read().splitlines()
    seq_dict = _read_FASTA_files(chr_paths, MOL_TYPE_CHR)
    pls_paths = open(pls_paths_file,'r').read().splitlines()
    seq_dict.update(_read_FASTA_files(pls_paths, MOL_TYPE_PLS))

    # Converting to dataframe
    seq_df = pd.DataFrame.from_dict(seq_dict).T
    seq_df.index.name='id'
    try:
        seq_df.to_csv(out_csv_file)
    except Exception as e:
        process_exception(f'Writing GC content  file {out_csv_file}: {e}')
    # Violinplot from dataframe
    second_df = pd.DataFrame(seq_df.to_dict())
    fig, ax = plt.subplots(figsize=(10,8))
    plt.suptitle('GC content ratio distribution')
    ax = sns.violinplot(x=MOL_TYPE_KEY, y=GC_RATIO_KEY, data=second_df)
    try:
        plt.savefig(out_png_file)
    except Exception as e:
        process_exception(f'Writing GC content  file {out_png_file}: {e}')
    # Determining equally spaced intervals from dataframe
    pls_df = seq_df[seq_df[MOL_TYPE_KEY] == MOL_TYPE_PLS]
    max_gc = pls_df[GC_RATIO_KEY].max()
    min_gc = pls_df[GC_RATIO_KEY].min()
    gc_bin_size = (max_gc - min_gc) / (n_gcints-2)
    gc_endpoints = (np.linspace(min_gc - gc_bin_size/2, max_gc + gc_bin_size/2, n_gcints - 1))
    try:
        with open(out_intervals_file, "w") as f:
            intervals_str = '\n'.join([str(x) for x in gc_endpoints])
            f.write(f'0\n{intervals_str}\n1\n')
    except Exception as e:
        process_exception(
            f'GC\tWriting {out_intervals_file}: {e}'
        )
            
def read_gc_intervals_file(gc_intervals_file):
    """
    Read GC intervals file
    Args:
        - gc_intervals_file (str): path to GC intervals file
    Returns:
        List(float): sorted list of GC intervals boundaries
    """
    with open(gc_intervals_file) as in_file:
        intervals = [
            float(x.rstrip())
            for x in in_file.readlines()
        ]
    
    check_number_eq(
        intervals[0], 0.0,
        msg=f'GC\t{gc_intervals_file}, First interval boundary must be 0, {intervals[0]} read'
    )
    check_number_eq(
        intervals[-1], 1.0,
        msg=f'GC\t{gc_intervals_file}, Last interval boundary must be 1, {intervals[-1]} read'
    )
    for x in range(1,len(intervals)):
        check_number_lt(
            intervals[x-1], intervals[x],
            msg=f'Intervals boundaries must be increasing {intervals[x]} <= {intervals[x-1]}'
        )
    return intervals

def intervals_boundaries_to_str(intervals):
    """
    Convert a list of intervals boundaries into a list of strings
    describing intervals of the form <start>-<end>
    """
    return [
        f'{intervals[i]}-{intervals[i+1]}'
        for i in range(len(intervals) - 1)
    ]

def compute_gc_probabilities(contigs_gc, gc_intervals):
    """
    Computes GC probabilities for a set of contigs and GC intervals

    Args:
        - contigs_gc (Dictionary): contig id -> {
          GC_COUNT_KEY: GC count, GC_RATIO_KEY: GC content ratio, LENGTH_KEY: contig length
          }
        - gc_intervals (List(float)): GC intervals boundaries

    Returns:
        (Dictionary) contig id -> List(float)
        where element in position i in the list is the probability to be in the (i+1)th interval
    """
    def _gprob2(n, g, p, m):
        """
        Compute probability of observing g GC nucleotides in a contig
        of length n within a molecule of GC content p using pseducount m.
        Done via logarithm to avoid overflow.
        """
        def _combln(n, k):
            """ Compute ln of n choose k. Note than n! = gamma(n+1) """
            return sc.gammaln(n + 1) - sc.gammaln(k + 1) - sc.gammaln(n - k + 1)
        alpha = m*p
        beta = m*(1-p)
        resultln = _combln(n, g) + sc.betaln(g + alpha, n - g + beta) - sc.betaln(alpha, beta)
        return math.exp(resultln)
    
    m = 10
    ctgs_gcp = {}
    for ctg_id,ctg_data in contigs_gc.items():
        n = ctg_data[LENGTH_KEY]
        g = ctg_data[GC_COUNT_KEY]
        total = 0
        gcp_array = []
        for i in range(0, len(gc_intervals)-1):
            gp2 = integrate.quad(
                lambda x: _gprob2(n,g,x,m),
                    gc_intervals[i], gc_intervals[i+1]
            )
            gp2 = gp2[0]/(gc_intervals[i+1] - gc_intervals[i])
            total += gp2
            gcp_array.append(gp2)
        ctgs_gcp[ctg_id] = [gcp/total for gcp in gcp_array]
    return ctgs_gcp

def _write_gc_probabilities_file(ctg_gcp_dict, gcp_out_file):
    try:
        with open(gcp_out_file, 'w') as out_file:
            for ctg_id,ctg_gcp in ctg_gcp_dict.items():
                gcp_str = '\t'.join([str(gcp) for gcp in ctg_gcp])
                out_file.write(f'{ctg_id}\t{gcp_str}\n')
    except Exception as e:
        process_exception(f'Writing GC probabilities file {gcp_out_file}: {e}')
                
def read_gc_probabilities_file(gcp_in_file, num_intervals=DEFAULT_N_GCINTS):
    """
    Args:
        - gcp_in_file (str): path to a GC probabilities file
        - num_intervals (int): expected number of GC content intervals
    Returns 
        Dictionary contig id -> list of probabilities ordered by GC interval
    """
    gcp_dict = {}
    with open(gcp_in_file, 'r') as in_file:
        for gcp_line in in_file.readlines():
            line = gcp_line.rstrip()
            ctg_data = line.split()
            check_num_fields_eq(
                ctg_data, num_intervals+1,
                msg=f'GC\t{gcp_in_file} {line}'
            )
            ctg_id = ctg_data[0]
            ctg_gcp_list = []
            for x in range(1,num_intervals+1):
                gcp = float(ctg_data[x])
                check_number_range(
                    gcp, (0.0,1.0),
                    msg=f'GC\t{gcp_in_file} {ctg_id} {gcp} probability'
                )
                ctg_gcp_list.append(gcp)
            ctg_gcp_total = round(sum(ctg_gcp_list),5)
            check_number_range(
                ctg_gcp_total, (0.0,1.0),
                msg=f'GC\t{gcp_in_file} {ctg_id} cumulated probability'
            )
            gcp_dict[ctg_id] = ctg_gcp_list
    return gcp_dict
            
def compute_gc_probabilities_file(gfa_file, gc_intervals_file, gcp_out_file, gfa_gzipped=True):
    """
    Computes probability that contigs originate from a molecule of a given GC content ratio,
    for a given list of GC content intervals.

    Args:
        - gfa_file (str): path to an  GFA file
        - gc_intervals_file (str) path to a GC content ratio intervals file
        - out_file: output file
        - gfa_gzipped (bool): True if input GFA file is gzipped

    Returns:
        creates out_file in format
        line 1: CTG<TAB><TAB-separated list of GC intervals>
        line 2+: GFA contig id<TAB><TAB-separated list of probabilities the contig oiginates from 
                 a molecule of GC content ratio within the interval>
    """
    gc_intervals = (
        DEFAULT_GC_INTERVALS if gc_intervals_file is None
        else read_gc_intervals_file(gc_intervals_file)
    )
    ctg_gc = read_GFA_ctgs(
        gfa_file,
        [GFA_SEQ_KEY],
        gzipped=gfa_gzipped,
        ctg_fun=_process_GFA_ctg
    )
    ctg_gc_probabilities = compute_gc_probabilities(ctg_gc, gc_intervals)
    _write_gc_probabilities_file(ctg_gc_probabilities, gcp_out_file)
   
