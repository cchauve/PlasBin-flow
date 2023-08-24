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

# Generic file functions

def __open_file_read(in_file_path, gzipped=False):
    """
    Open a file for reading

    Args:
        - in_file_path (str): path of file to read
        - gzipped (bool): True if file is gzipped

    Returns:
        - (file object) file to read
    """
    if gzipped:
        return gzip.open(in_file_path, 'rt')
    else:
        return open(in_file_path, 'r')

def __open_file_write(in_file_path, gzipped=False):
    """
    Open a file for writing

    Args:
        - in_file_path (str): path of file to write to
        - gzipped (bool): True if file is gzipped

    Returns:
        - (file object) file to write to
    """
    if gzipped:
        return gzip.open(in_file_path, 'wt')
    else:
        return open(in_file_path, 'w')
    
## Reading FASTA files

def gunzip_FASTA(in_file_path, out_file_path):
    """
    Gunzip a FASTA file

    Args:
       in_file_path (str): path to input gzipped FASTA file
       out_file_path (str): path to output FASTA file

    Returns:
       Creates FASTA file out_file_path
    """
    records = []
    with gzip.open(in_file_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
        with open(out_file_path, 'w') as out_file:
            SeqIO.write(records, out_file, 'fasta')

def read_FASTA_ctgs(in_file_path, record_fun, gzipped=False):
    """
    Read FASTA file, processing each sequence

    Args:
        - in_file_path (str): path of FASTA file to read
        - record_fun (function): processing function taking a single input of type SeqIO
        - gzipped (bool): True if file is gzipped

    Returns:
        - (Dictionary) sequence id (str) -> record_fun(sequence.record)
    """
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

def read_FASTA_id(in_file_path, gzipped=False):
    """
    Computes the list of sequences id in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - gzipped (bool): True if file is gzipped

    Returns:
        - (List(str)): list of sequence ids
    """
    return list(
        read_FASTA_ctgs(
            in_file_path,
            record_fun=lambda x: None,
            gzipped=gzipped
        ).keys()
    )

def read_FASTA_len(in_file_path, gzipped=False):
    """
    Computes the length of entry sequences in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - gzipped (bool): True if file is gzipped

    Returns:
        - (Dictionary): sequence id (str) -> length of sequence (int)
    """
    return read_FASTA_ctgs(
        in_file_path,
        record_fun=lambda x: len(x.seq),
        gzipped=gzipped
    )

def read_FASTA_seq(in_file_path, gzipped=False):
    """
    Computes entry sequences in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - gzipped (bool): True if file is gzipped

    Returns:
        - (Dictionary): sequence id (str) -> sequence (str)
    """
    return read_FASTA_ctgs(
        in_file_path,
        record_fun=lambda x: x.seq,
        gzipped=gzipped
    )

## Reading GFA files (contigs and links only)

def gunzip_GFA(in_file_path, out_file_path):
    """
    Gunzip a GFA file

    Args:
       in_file_path (str): path to input gzipped GFA file
       out_file_path (str): path to output GFA file

    Returns:
       Creates oGFA file out_file_path
    """
    with gzip.open(in_file_path) as in_file, open(out_file_path, 'wb') as out_file:
        shutil.copyfileobj(in_file, out_file)

# Mandatory fields in GFA contigs and links
GFA_SEQ_KEY = 'Sequence'
GFA_FROM_KEY = 'From'
GFA_FROM_ORIENT_KEY = 'FromOrient'
GFA_TO_KEY = 'To'
GFA_TO_ORIENT_KEY = 'ToOrient'
GFA_OVERLAP_KEY = 'Overlap'

# Attributes types as defined in http://gfa-spec.github.io/GFA-spec/GFA1.html
GFA_DEFAULT_TYPES = {
    GFA_SEQ_KEY: lambda x: str(x),
    'LN': lambda x: int(x),
    'RC': lambda x: int(x),
    'FC': lambda x: int(x),
    'KC': lambda x: int(x),
    'SH': lambda x: str(x),
    'UR': lambda x: str(x),
    GFA_FROM_KEY: lambda x: str(x),
    GFA_FROM_ORIENT_KEY: lambda x: str(x),
    GFA_TO_KEY: lambda x: str(x),
    GFA_TO_ORIENT_KEY: lambda x: str(x),
    GFA_OVERLAP_KEY: lambda x: str(x),
    'MQ': lambda x: int(x),
    'MN': lambda x: int(x),
    'ID': lambda x: str(x)
}
        
def __add_attributes(attributes_data, attributes_fun):
    """
    Creates a dictionary of attributes for a contig/link

    Args:
        - attributes_data (List): list of attributes in format str(key:type:value)
        - attributes_fun (Dictionary): attribute key (str) -> conversion/processing function

    Returns:
       - (Dictionary) attribute key: attribute value (None if missing attribute)
    """
    attributes_list = list(attributes_fun.keys())
    attributes_dict = {att_key: None for att_key in attributes_list}
    for att_data in attributes_data:
        att_key = att_data.split(':')[0]
        att_val = att_data.split(':')[2]
        if att_key in attributes_list:
            attributes_dict[att_key] = attributes_fun[att_key](att_val)
    return attributes_dict

def read_GFA_ctgs(in_file_path, attributes_fun, gzipped=False):
    """
    Read contigs and their attributes from a GFA files

    Args:
        - in_file_path (str): path to GFA file to read
        - attributes_fun (Dictionary): attribute key (str) -> conversion/processing function
        - gzipped (bool): True if gzipped GFA file

    Returns:
       - (Dictionary) contig id -> (Dictionary) attribute key: attribute value (None if missing attribute)
         where attribute key GFA_SEQ_KEY is for the contig sequence
    Assumption: 
       - every contig has an associated id and sequence (not checked)
    """
    result = {}
    with __open_file_read(in_file_path, gzipped) as in_file:
        for line in [x for x in in_file.readlines() if x[0]=='S']:
            ctg_data = line.strip().split('\t')
            ctg_id,ctg_seq = ctg_data[1],ctg_data[2]
            result[ctg_id] = __add_attributes(
                [f'{GFA_SEQ_KEY}::{ctg_seq}'] + ctg_data[3:],
                attributes_fun
            )
    return result

def read_GFA_links(in_file_path, attributes_fun, gzipped=False):
    """
    Read links and their attributes from a GFA files

    Args:
        - in_file_path (str): path to GFA file to read
        - attributes_fun (Dictionary): attribute key (str) -> conversion/processing function
        - gzipped (bool): True if gzipped GFA file

    Returns:
       - (Dictionary) contig id -> 
         List(links from contig id
           (Dictionary) attribute key: attribute value (None if missing attribute))
           graph attributes: GFA_FROM_ORIENT_KEY, GFA_TO_KEY, GFA_TO_ORIENT_KEY, GFA_OVERLAP_KEY

    Assumption: 
       - every contig has an associated id and sequence (not checked)
    """
    result = defaultidct(list)
    with __open_file_read(in_file_path, gzipped) as in_file:
        for line in [x for x in in_file.readlines() if x[0]=='L']:
            ctg_data = line.strip().split('\t')
            ctg_from = ctg_data[1]
            ctg_from_orient = ctg_data[2]
            ctg_to = ctg_data[3]
            ctg_to_orient = ctg_data[4]
            overlap = ctg_data[5]
            result[ctg_from].append(
                __add_attributes([
                    f'{GFA_TO_KEY}::{ctg_to}',
                    f'{GFA_FROM_ORIENT_KEY}::{ctg_from_orient}',
                    f'{GFA_TO_ORIENT_KEY}::{ctg_to_orient}',
                    f'{GFA_OVERLAP_KEY}::{overlap}'
                ] + ctg_data[6:], attributes_fun)
            )
    return result

def write_GFA_to_FASTA(in_GFA_file, out_FASTA_file, in_gzipped, out_gzipped):
    GFA_ctg_seqs = read_GFA_ctgs(
        in_GFA_file,
        attributes_fun={GFA_SEQ_KEY: GFA_DEFAULT_TYPES[GFA_SEQ_KEY]},
        gzipped=in_gzipped
    )
    ctg_records = [
        SeqRecord(Seq(y[GFA_SEQ_KEY]), id=x, name=x, description=f'{x}.GFA')
        for x,y in GFA_ctg_seqs.items()
    ]
    with __open_file_write(out_FASTA_file, gzipped=out_gzipped) as out_file:
        SeqIO.write(ctg_records, out_file, 'fasta')


