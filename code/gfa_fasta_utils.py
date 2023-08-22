#!/usr/bin/env python

"""
Manipulating FASTA and GFA files
"""

import csv
import gzip
import os
import shutil
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Reading FASTA files
def __open_file_read(in_file_path, gzipped=False):
    """
    Open in_file_path for reading
    """
    if gzipped:
        return gzip.open(in_file_path, 'rt')
    else:
        return open(in_file_path, 'r')

def __open_file_write(in_file_path, gzipped=False):
    """
    Open in_file_path for writing
    """
    if gzipped:
        return gzip.open(in_file_path, 'wt')
    else:
        return open(in_file_path, 'w')
    
def __read_FASTA(in_file_path, record_fun, gzipped=False):
    """
    Read FASTA file in_file_path
    Process each contig using record_fun
    Returns Dictionary sequence id -> record_fun(sequence.record)
    """
    if os.path.isfile(in_file_path):
        return {
            id: record_fun(record)
            for id,record in SeqIO.to_dict(
                    SeqIO.parse(
                        __open_file_read(
                            in_file_path,
                            gzipped
                        ),
                        'fasta'
                    )
            ).items()
        }
    else:
        return {}

def read_FASTA_id(in_file_path, gzipped=False):
    """
    Computes the list of sequences id
    """
    return list(
        __read_FASTA(
            in_file_path,
            record_fun=lambda x: None,
            gzipped=gzipped
        ).keys()
    )

def read_FASTA_len(in_file_path, gzipped=False):
    """
    Computes Dictionary sequence id -> sequence length
    """
    return __read_FASTA(
        in_file_path,
        record_fun=lambda x: len(x.seq),
        gzipped=gzipped
    )

def read_FASTA_seq(in_file_path, gzipped=False):
    """
    Computes Dictionary sequence id -> sequence 
    """
    return __read_FASTA(
        in_file_path,
        record_fun=lambda x: x.seq,
        gzipped=gzipped
    )

## Reading GFA files
## TODO: read edges

def __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: x, gzipped=False):
    """
    Read a GFA files and process its contigs
    Reaturns Dictionary contig id -> ctg_fun(contig data)
    """
    result = {}
    if os.path.isfile(in_file_path):
        with __open_file_read(in_file_path, gzipped) as in_file:
            for line in in_file.readlines():
                line_split = line.strip().split('\t')
                if line[0] == 'S':
                    ctg_id,ctg_data = line_split[1],line_split[2:]
                    result[ctg_id] = ctg_fun(ctg_data)
    return result

def __GFA_S_attributes(in_data, in_key):
    attributes = {
        x.split(':')[0]: x.split(':')[2]
        for x in in_data
    }
    if in_key in attributes.keys():
        return attributes[in_key]
    else:
        return None

def read_GFA_ctg_id(in_file_path, gzipped=False):
    """
    Computes the list of contigs id
    """
    return list(
        __read_GFA_ctgs(
            in_file_path,
            ctg_fun=lambda x: None,
            gzipped=gzipped
        ).keys()
    )
    
def read_GFA_ctg_len(in_file_path, gzipped=False):
    """
    Computes Dictionary contig id -> contig length
    """
    return __read_GFA_ctgs(
        in_file_path,
        ctg_fun=lambda x: len(x[0]),
        gzipped=gzipped
    )

def read_GFA_ctg_seq(in_file_path, gzipped=False):
    """
    Computes Dictionary contig id -> contig sequence
    """
    return __read_GFA_ctgs(
        in_file_path,
        ctg_fun=lambda x: x[0],
        gzipped=gzipped
    )

def read_GFA_ctg_RC(in_file_path, gzipped=False):
    """
    Computes Dictionary contig id -> contig RC value (None if absent)
    """    
    return __read_GFA_ctgs(
        in_file_path,
        ctg_fun=lambda x: int(__GFA_S_attributes(x, 'RC')),
        gzipped=gzipped
    )

def read_GFA_ctg_KC(in_file_path, gzipped=False):
    """
    Computes Dictionary contig id -> contig KC value (None if absent)
    """    
    return __read_GFA_ctgs(
        in_file_path,
        ctg_fun=lambda x: int(__GFA_S_attributes(x, 'KC')),
        gzipped=gzipped
    )

def read_GFA_ctg_LN(in_file_path, gzipped=False):
    """
    Computes Dictionary contig id -> contig LN value (None if absent)
    """    
    return __read_GFA_ctgs(
        in_file_path,
        ctg_fun=lambda x: int(__GFA_S_attributes(x, 'LN')),
        gzipped=gzipped
    )

def read_GFA_ctg_dp(in_file_path, gzipped=False):
    """
    Computes Dictionary contig id -> contig dp value (None if absent)
    """    
    return __read_GFA_ctgs(
        in_file_path,
        ctg_fun=lambda x: float(__GFA_S_attributes(x, 'dp')),
        gzipped=gzipped
    )

def write_GFA_to_FASTA(in_GFA_file, out_FASTA_file, in_gzipped, out_gzipped):
    GFA_ctg_seqs = read_GFA_ctg_seq(in_GFA_file, gzipped=in_gzipped)
    ctg_records = [
        SeqRecord(Seq(y), id=x, name=x, description=f'{x}.GFA')
        for x ,y in GFA_ctg_seqs.items()
    ]
    with __open_file_write(out_FASTA_file, gzipped=out_gzipped) as out_file:
        SeqIO.write(ctg_records, out_file, 'fasta')

def gunzip_FASTA(in_fasta_file, out_fasta_file):
    """
    Gunzip a FASTA file

    Args:
       in_fasta_file (str): path to input gzipped FASTA file
       out_fasta_file (str): path to output FASTA file

    Returns:
       Creates out_fasta_file or log a warning if it already exists
    """
    records = []
    with gzip.open(in_fasta_file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
        with open(out_fasta_file, 'w') as out_file:
            SeqIO.write(records, out_file, 'fasta')

def gunzip_GFA(in_gfa_file, out_gfa_file):
    """
    Gunzip a GFA file

    Args:
       in_gfa_file (str): path to input gzipped GFA file
       out_gfa_file (str): path to output GFA file

    Returns:
       Creates out_gfa_file or log a warning if it already exists
    """
    with gzip.open(in_gfa_file) as in_file, open(out_gfa_file, 'wb') as out_file:
        shutil.copyfileobj(in_file, out_file)
