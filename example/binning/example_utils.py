""" Tools to create input and parse output for classification external tools """

import sys
import os
sys.path.append('../../code/')
from collections import defaultdict
import csv

from Bio import SeqIO

from gfa_fasta_utils import (
    write_GFA_to_FASTA,
    read_FASTA_len
)

def create_FASTA_file(in_GFA_file, out_FASTA_file):
    """
    Create a FASTA file from a gzipped GFA file

    Args:
        - in_GFA_file (str): path to gzipped GFA file
        - out_FASTA_file (str): path to unzipped FASTA file

    Returns: Creates file out_FASTA_file
    """
    write_GFA_to_FASTA(in_GFA_file, out_FASTA_file, in_gzipped=True, out_gzipped=False)

def filter_FASTA_file(in_FASTA_file, out_FASTA_file, min_len):
    """
    Filter contigs from input FASTA to discard contigs of length
    below min_len and write the resulting file in out_FASTA_file
    """
    with open(in_FASTA_file) as in_file:
        records = []
        for record in SeqIO.parse(in_file, 'fasta') :
            if len(record.seq) >= min_len :
                records.append(record)
    with open(out_FASTA_file, 'w') as out_file:
        SeqIO.write(records, out_file, 'fasta')

def _read_RFPlasmid_results(in_RFPlasmid_file):
    """
    Read sample RFPlasmid results file in_RFPlasmid_file
    Returns dictionary: contig -> plasmid score
    """
    sample_scores_dict = {}
    with open(in_RFPlasmid_file) as csvfile:
        csv_reader = csv.DictReader(csvfile, delimiter=',')
        for ctg_result in csv_reader:
            ctg_id = ctg_result['contigID'].split()[0]
            pls_score = float(ctg_result['votes plasmid'])
            pred = ctg_result['prediction']
            sample_scores_dict[ctg_id] = pls_score
    return sample_scores_dict

def _read_PlasClass_results(in_PlasClass_file):
    """
    Read sample PlasClass results file in_PlasClass_file
    Returns dictionary: contig -> plasmid score
    """
    sample_scores_dict = {}
    with open(in_PlasClass_file) as csvfile:
        csv_reader = csv.DictReader(
            csvfile, fieldnames=['contig', 'score'], delimiter='\t'
        )
        for ctg_result in csv_reader:
            ctg_id = ctg_result['contig']
            pls_score = float(ctg_result['score'])
            sample_scores_dict[ctg_id] = pls_score
    return sample_scores_dict

def _read_plASgraph_results(in_plASgraph_file):
    """
    Read sample plASgraph results file in_plASgraph_file
    Returns dictionary: contig -> plasmid score
for contigs of sample sample_id
    """
    sample_scores_dict = {}
    with open(in_plASgraph_file) as csvfile:
        csv_reader = csv.DictReader(csvfile, delimiter=',')
        for ctg_result in csv_reader:
            sample = ctg_result['sample']
            ctg_id = ctg_result['contig']
            pls_score = float(ctg_result['plasmid_score'])
            pred = ctg_result['label'][0]
            sample_scores_dict[ctg_id] = pls_score
    return sample_scores_dict

def reformat_tool_output(in_tool_file, read_tool_fun, tool_name, out_scores_file):
    """
    Reformat an RFPlasmid file for one sample

    Args:
        - in_tool_file (str): path to the result file of external tool
        - read_tool_fun (function): read the external tool results
        - tool_name (str): external tool name
        - out_scores_file (str): path to created file

    Returns: Creates files out_file
        format:
        contig_id<TAB>score
    """
    sample_scores_dict = read_tool_fun(in_tool_file)
    with open(out_scores_file, 'w') as out_file:
        for ctg_id,ctg_data in sample_scores_dict.items():
            score = sample_scores_dict[ctg_id]
            out_file.write(f'{ctg_id}\t{score}\n')
        
def main():
    cmd = sys.argv[1]

    if cmd == 'create_FASTA':
        in_GFA_file = sys.argv[2]
        out_FASTA_file = sys.argv[3]
        create_FASTA_file(in_GFA_file, out_FASTA_file)

    elif cmd == 'filter_FASTA':
        in_FASTA_file = sys.argv[2]
        out_FASTA_file = sys.argv[3]
        min_len = int(sys.argv[4])
        filter_FASTA_file(in_FASTA_file, out_FASTA_file, min_len)
    
    elif cmd == 'reformat_RFPlasmid':
        in_RFPlasmid_file = sys.argv[2]
        out_file = sys.argv[3]
        reformat_tool_output(
            in_RFPlasmid_file, _read_RFPlasmid_results, 'RFPlasmid', out_file
        )

    elif cmd == 'reformat_PlasClass':
        in_PlasClass_file = sys.argv[2]
        out_file = sys.argv[3]
        reformat_tool_output(
            in_PlasClass_file, _read_PlasClass_results, 'PlasClass', out_file
        )

    elif cmd == 'reformat_plASgraph':
        in_plASgraph_file = sys.argv[2]
        out_file = sys.argv[3]
        reformat_tool_output(
            in_plASgraph_file, _read_plASgraph_results, 'plASgraph', out_file
        )

if __name__ == "__main__":
    main()
