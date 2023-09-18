""" Functions to compute the gene density of contigs """

from gfa_fasta_utils import (
    read_GFA_len
)
from mappings_utils import (
    read_blast_outfmt6_file,
    compute_blast_s_intervals,
    filter_blast_outfmt6
)
from log_errors_utils import (
    process_exception,
    check_number_range,
    check_num_fields
)

GD_INTERVALS_KEY = 'intervals'
GD_DENSITY_KEY = 'gd'

def compute_gene_density(
    gfa_file, mappings_file,
    pid_threshold, cov_threshold,
    gfa_gzipped=True
):
    """
    Computes gene density and covering intervals for all contigs

    Args:
        - gfa_file (str): path to a GFA file
        - mappings_file (str): path to a genes to contigs mappings file
        - pid_threshold (float): min identity percentage to kep a mapping
        - cov_threshold (float): min gene coverage to keep a mapping
        - gfa_gzipped (bool): True if GFA file is gzipped

    Returns:
        (Dictionary): contig id (str) ->
                      (Dictionary): GD_DENSITY_KEY: gene density (float)
                                    GD_INTERVALS_KEY: List((gene,start,end)) list of covering intervals
    """
    def _get_union(intervals):
        """
        Takes the gene covering intervals for a contig and finds their union
        The length of the union is used to compute gene coverage
        """
        intervals_union = []
        for _,start,end in intervals:
            if intervals_union and intervals_union[-1][1] >= start-1:
                intervals_union[-1][1] = max(intervals_union[-1][1],end)
            else:
                intervals_union.append([start,end])
        return intervals_union

    def _compute_gd(intervals_union, ctg_len):
        """
        Computes gene density using list of coverage intervals and contig length
        """
        covered = 0
        for interval in intervals_union:
            covered += interval[1] - interval[0] + 1
        return covered / ctg_len
    
    mappings_df = read_blast_outfmt6_file(mappings_file)
    filter_blast_outfmt6(mappings_df, min_pident=pid_threshold, min_q_cov=cov_threshold)
    ctg_len = read_GFA_len(gfa_file, gzipped=gfa_gzipped)
    ctg_intervals = compute_blast_s_intervals(mappings_df)
    ctg_gd_dict = {
        ctg_id: {GD_DENSITY_KEY: 0.0, GD_INTERVALS_KEY: []}
        for ctg_id in ctg_len.keys()
    }
    for ctg_id,intervals in ctg_intervals.items():
        sorted_intervals = sorted(intervals, key=lambda x: x[1])
        intervals_union = _get_union(sorted_intervals)
        ctg_gd = _compute_gd(
            intervals_union, ctg_len[ctg_id]
        )
        ctg_gd_dict[ctg_id] = {
            GD_DENSITY_KEY: ctg_gd,
            GD_INTERVALS_KEY: intervals_union
        }
    return ctg_gd_dict

def _write_intervals(intervals_list):
    return ' '.join([
        f'{start}:{end}'
        for [start,end] in intervals_list
    ])

def _read_intervals(intervals_str):
    return [
        [
            int(interval.split(':')[1]),
            int(interval.split(':')[2])
        ]
        for interval in intervals_str.split()
    ]

def _write_gene_density_file(ctg_gd_dict, gd_out_file):
    """
    Writes the gene density of contigs of a sample in a TSV file

    Args:
        - ctg_gd_dict (Dictionary): gene density dictionary (see compute_gene_density)
        - gd_out_file (str): path to the file to be written

    Returns:
        Creates gd_out_file with format
        <contig id><TAB><gene density><TAB><space separated intervals <gene id>:<start>:<end>
    """
    try:
        with open(gd_out_file, 'w') as out_file:
            for ctg_id,ctg_data in ctg_gd_dict.items():
                ctg_gd = ctg_data[GD_DENSITY_KEY]
                ctg_intervals_str = _write_intervals(ctg_data[GD_INTERVALS_KEY])
                out_file.write(f'{ctg_id}\t{ctg_gd}\t{ctg_intervals_str}\n')
    except Exception as e:
        process_exception(f'Writing gene density file {gd_out_file}: {e}')

def compute_gene_density_file(
        gfa_file, mappings_file, gd_out_file,
        pid_threshold, cov_threshold,
        gfa_gzipped=True
):
    """
    Computes gene density and covering intervals for all contigs, writes in a TSV file

    Args:
        - gfa_file (str): path to a GFA file
        - mappings_file (str): path to a genes to contigs mappings file
        - gd_out_file (str): file where to write gene density
        - pid_threshold (float): min identity percentage to kep a mapping
        - cov_threshold (float): min gene coverage to keep a mapping
        - gfa_gzipped (bool): True if GFA file is gzipped

    Returns:
        Creates gd_out_file, format see _write_gene_density_file    
    """
    ctg_gd_dict = compute_gene_density(
        gfa_file, mappings_file,
        pid_threshold, cov_threshold,
        gfa_gzipped=gfa_gzipped
    )
    _write_gene_density_file(ctg_gd_dict, gd_out_file)

def read_gene_density_file(gd_file, read_intervals=False):
    """
    Reads a gene density file

    Args:
        - gd_file (str): path to the gene density file in format 
         <contig id><TAB><gene density>[OPTIONAL<TAB><space separated intervals <gene id>:<start>:<end>>]
        - read_intervals (bool): if False, does not read covering intervals

    Returns:
        - if read_intervals: (Dictionary) as in compute_gene_density
        - else: (Dictionary) contig id -> gene density (float)
    """
    ctg_gd_dict = {}
    with open(gd_file, 'r') as in_file:
        for ctg_line in in_file.readlines():
            line = ctg_line.rstrip()
            ctg_data = line.rstrip().split('\t')
            check_num_fields(
                ctg_data, 2,
                msg=f'GD\t{gd_file} {line}'
            )
            ctg_id = ctg_data[0]
            ctg_gd = float(ctg_data[1])
            check_number_range(
                ctg_gd, (0.0,1.0),
                msg=f'GD\t{gd_file} {line}'
            )
            if read_intervals:
                check_num_fields(
                    ctg_data, 3,
                    msg=f'GD\t{gd_file} {line}'
                )
                ctg_intervals = _read_intervals(ctg_data[2])
                ctg_gd_dict[ctg_id] = {
                    GD_DENSITY_KEY: ctg_gd,
                    GD_INTERVALS_KEY: ctg_intervals
                }
            else:
                ctg_gd_dict[ctg_id] = ctg_gd
    return ctg_gd_dict
