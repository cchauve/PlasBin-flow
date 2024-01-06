"""
Manipulating FASTA and GFA files
"""

import gzip
import os
import shutil
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from log_errors_utils import (
    process_exception,
    check_num_fields,
    CustomException
)

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

def read_FASTA_ctgs(in_file_path, record_fun, gzipped=False, id_fun=lambda x: x):
    """
    Read FASTA file, processing each entry

    Args:
        - in_file_path (str): path of FASTA file to read
        - record_fun (function): processing function taking a single input of type SeqIO
        - gzipped (bool): True if file is gzipped
        - id_fun (function) processing function for record is

    Returns:
        - (Dictionary) sequence id (str) -> record_fun(sequence.record)
    """
    try:
        ctgs_dict = {
            id_fun(id): record_fun(record)
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
    except Exception as e:
        process_exception(f'Reading {in_file_path}: {e}')
    return ctgs_dict
        

def read_FASTA_id(in_file_path, gzipped=False, id_fun=lambda x: x):
    """
    Computes the list of sequences id in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - gzipped (bool): True if file is gzipped
        - id_fun (function) processing function for record is

    Returns:
        - (List(str)): list of sequence ids
    """
    return [
        id_fun(ctg_id)
        for ctg_id in list(
                read_FASTA_ctgs(
                    in_file_path,
                    record_fun=lambda x: None,
                    gzipped=gzipped
                ).keys()
        )
    ]

def read_FASTA_len(in_file_path, gzipped=False, id_fun=lambda x: x):
    """
    Computes the length of entry sequences in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - gzipped (bool): True if file is gzipped
        - id_fun (function) processing function for record is

    Returns:
        - (Dictionary): sequence id (str) -> length of sequence (int)
    """
    return read_FASTA_ctgs(
        in_file_path,
        record_fun=lambda x: len(x.seq),
        gzipped=gzipped,
        id_fun=id_fun
    )

def read_FASTA_seq(in_file_path, gzipped=False, id_fun=lambda x: x):
    """
    Computes entry sequences in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - gzipped (bool): True if file is gzipped
        - id_fun (function) processing function for record is

    Returns:
        - (Dictionary): sequence id (str) -> sequence (str)
    """
    return read_FASTA_ctgs(
        in_file_path,
        record_fun=lambda x: str(x.seq),
        gzipped=gzipped,
        id_fun=id_fun
    )

## Reading GFA files (contigs and links only)

def gunzip_GFA(in_file_path, out_file_path):
    """
    Gunzip a GFA file

    Args:
       in_file_path (str): path to input gzipped GFA file
       out_file_path (str): path to output GFA file

    Returns:
       Creates GFA file out_file_path
    """
    try:
        with gzip.open(in_file_path) as in_file, open(out_file_path, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)
    except Exception as e:
        process_exception(f'FASTA\tGunzipping {in_file_path} to {out_file_path}: {e}')
            
# Mandatory fields in GFA contigs and links
GFA_SEQ_KEY = 'Sequence'
GFA_LEN_KEY = 'Length'
GFA_FROM_KEY = 'From'
GFA_FROM_ORIENT_KEY = 'FromOrient'
GFA_TO_KEY = 'To'
GFA_TO_ORIENT_KEY = 'ToOrient'
GFA_OVERLAP_KEY = 'Overlap'

# Conversionto of GFA attributes.
# Missing attributes types: B, J
GFA_ATTRIBUTE_TYPE = {
    'i': lambda x: int(float(x)),
    'f': lambda x: float(x),
    'Z': lambda x: str(x),
    'A': lambda x: str(x),
    'H': lambda x: bytes(x)
}
# 'i' converted into float at first to handle cases such as '1.20934e+06'

def __add_attributes(attributes_data, attributes_list):
    """
    Creates a dictionary of attributes for a contig/link

    Args:
        - attributes_data (List): list of attributes in format str(key:type:value)
        - attributes_list (List(str)): list of attribute keys to read
          ['all'] for recording all attributes 

    Returns:
        - (Dictionary) attribute key: attribute value (None if missing attribute)
        attributes values are converted to the type defined by GFA_ATTRIBUTE_TYPE
        if an attribute key in attributes_list is not a key in GFA_ATTRIBUTE_TYPE
        the attribute value is recorded as a string
    """
    attributes_dict = {
        att_key: None
        for att_key in attributes_list
        if attributes_list != ['all']
    }
    for att_data in attributes_data:
        att_split = att_data.split(':')
        check_num_fields(att_split, 3)
        att_key,att_type = att_split[0:2]
        att_val = ':'.join(att_split[2:])
        if attributes_list == ['all'] or att_key in attributes_list:
            if att_type not in GFA_ATTRIBUTE_TYPE.keys():
                attributes_dict[att_key] = att_val
            else:
                attributes_dict[att_key] = GFA_ATTRIBUTE_TYPE[att_type](att_val)
    return attributes_dict

def __write_attributes(attributes_dict, keys_to_remove=[], sep=' '):
    """
    Write GFA attributes into a string

    Args:
        - attributes_dict (Dictionary): attribute key -> attribute value
        - keys_to_remove (List(str)): list of keys to not print
        - sep (str): separating string

    Returns:
        (str): list of attributes in format sep.join(key:value)
        attributes with None value are not written
    """
    return sep.join(
        [
            f'{x}:{y}'
            for x,y in attributes_dict.items()
            if x not in keys_to_remove and y is not None
        ]
    )

def __assert_attributes_list(attributes_list):
    """
    Assert that an attributes list is either ['all'] or does not contain 'all'
    Used only for development

    Args:
        - attributes_list (List(str))
    """
    assert (
        attributes_list == ['all'] or 'all' not in attributes_list
    ), f'incorrect GFA attributes list {attributes_list}'

def read_GFA_ctgs(in_file_path, attributes_list, gzipped=False, ctg_fun=lambda x: x, id_fun=lambda x: x):
    """
    Read contigs and their attributes from a GFA files

    Args:
        - in_file_path (str): path to GFA file to read
        - attributes_list (List(str)): list of attribute keys to read
          ['all'] for recording all attributes
        - gzipped (bool): True if gzipped GFA file
        - ctg_fun: function that process a contig information
        - id_fun: function that process a contig id

    Returns:
       - (Dictionary) contig id -> ctg_fun(
           (Dictionary) attribute key: attribute value 
           (None if missing attribute)
         )
         where attribute key GFA_SEQ_KEY is for the contig sequence
    Assumption: 
       - every contig has an associated id and sequence (not checked)
    """
    __assert_attributes_list(attributes_list)
    result = {}
    with __open_file_read(in_file_path, gzipped) as in_file:
        for gfa_line in [x for x in in_file.readlines() if x[0]=='S']:
            line = gfa_line.rstrip()
            ctg_data = line.split('\t')
            check_num_fields(ctg_data, 2)
            ctg_id,ctg_seq = ctg_data[1],ctg_data[2]
            ctg_len = len(ctg_seq)
            att_data = [
                f'{GFA_SEQ_KEY}:Z:{ctg_seq}', f'{GFA_LEN_KEY}:i:{ctg_len}'
            ] + ctg_data[3:]
            result[id_fun(ctg_id)] = ctg_fun(
                __add_attributes(att_data, attributes_list)
            )
    return result

def read_GFA_id(in_file_path, gzipped=False, id_fun=lambda x: x):
    """
    Computes the list of segments (contigs) id in a GFA file

    Args:
        - in_file_path (str): path of GFA file to read
        - gzipped (bool): True if file is gzipped
        - id_fun: function that process a contig id

    Returns:
        - (List(str)): list of sequence ids
    """
    return [
        id_fun(ctg_id)
        for ctg_id in list(
                read_GFA_ctgs(
                    in_file_path,
                    attributes_list=[],
                    ctg_fun=lambda x: None,
                    gzipped=gzipped
                ).keys()
        )
    ]

def read_GFA_attribute(in_file_path, att_key, gzipped=False, id_fun=lambda x: x):
    """
    Computes the length of segments (contigs) in a GFA file

    Args:
        - in_file_path (str): path of GFA file to read
        - att_key (str): attribute key
        - gzipped (bool): True if file is gzipped
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> attribute value
    """
    return  {
        ctg_id: ctg_attributes[att_key]
        for ctg_id,ctg_attributes in read_GFA_ctgs(
                in_file_path, 
                [att_key],
                gzipped=gzipped,
                id_fun=id_fun
        ).items()
    }

def read_GFA_len(in_file_path, gzipped=False, id_fun=lambda x: x):
    """
    Computes the length of segments (contigs) in a GFA file

    Args:
        - in_file_path (str): path of GFA file to read
        - gzipped (bool): True if file is gzipped
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> length of sequence (int)
    """
    return read_GFA_attribute(
        in_file_path, GFA_LEN_KEY, gzipped=gzipped, id_fun=id_fun
    )

def read_GFA_seq(in_file_path, gzipped=False, id_fun=lambda x: x):
    """
    Computes segments (contigs) sequences in a GFA file

    Args:
        - in_file_path (str): path of GFA file to read
        - gzipped (bool): True if file is gzipped
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> sequence (str)
    """
    return read_GFA_attribute(
        in_file_path, GFA_SEQ_KEY, gzipped=gzipped, id_fun=id_fun
    )

def _ctgs_normalized_coverage(attributes_dict):
    """
    Computes the normalized coverage for a set of contigs
    Args:
        - attributes_dict (Dictionary): contig id -> attributes dictionary, including LN and KC keys
    Returnd:
        (Dictionary): contig id (str) -> normalized coverage (float)
    """
    # Total assembly coverage divided by k (k-mer value)
    total_coverage = sum([
        ctg_data['KC'] for _,ctg_data in attributes_dict.items()
    ])
    # Total assembly length
    total_length = sum([
        ctg_data[GFA_LEN_KEY] for _,ctg_data in attributes_dict.items()
    ])
    # Average assembly coverage: k * total_coverage / total_length
    # Contig total coverage = k * ctg.KC
    # Contig average coverage = k * ctg.KC / ctg.Length
    # Contig normalized coverage = (k * ctg.KC / ctg.Length) / (k * total_coverage / total_length)
    # = ctg.KC * total_length / ctg.LN * total_coverage
    return {
        ctg_id: (ctg_data['KC'] * total_length) / (ctg_data[GFA_LEN_KEY] * total_coverage)
        for ctg_id, ctg_data in attributes_dict.items()
    }

def read_GFA_normalized_coverage(in_file_path, cov_key=None, gzipped=False, id_fun=lambda x: x):
    """
    Computes normalized coverage of segments (contigs) in a GFA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - cov_key (str or None): attribute key that records normalized coverage
           if None, coverage is based on attributes KC and LN
        - gzipped (bool): True if file is gzipped
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> normalized coverage (float)
    """
    
    if cov_key is None:
        ctgs_data = read_GFA_ctgs(
            in_file_path,
            [GFA_LEN_KEY, 'KC'],
            gzipped=gzipped,
            id_fun=id_fun
        )
        try:
            result = _ctgs_normalized_coverage(ctgs_data)
        except Exception as e:
            process_exception(f'GFA\tComputing normalized coverage {in_file_path}: {e}')
        else:
            return result
    else:
        return read_GFA_attribute(
            in_file_path, cov_key, gzipped=gzipped, id_fun=id_fun
        )

def read_GFA_links(in_file_path, gzipped=False):
    """
    Read links and their attributes from a GFA files

    Args:
        - in_file_path (str): path to GFA file to read
        - gzipped (bool): True if gzipped GFA file

    Returns:
       - (Dictionary) contig id -> 
         List(links from contig id
           (Dictionary) attribute key: attribute value
           graph attributes: GFA_FROM_ORIENT_KEY, GFA_TO_KEY, GFA_TO_ORIENT_KEY, GFA_OVERLAP_KEY
    """
    result = defaultdict(list)
    with __open_file_read(in_file_path, gzipped) as in_file:
        for gfa_line in [x for x in in_file.readlines() if x[0]=='L']:
            line = gfa_line.rstrip()
            ctg_data = line.split('\t')
            check_num_fields(ctg_data, 5)
            ctg_from = ctg_data[1]
            ctg_from_orient = ctg_data[2]
            ctg_to = ctg_data[3]
            ctg_to_orient = ctg_data[4]
            overlap = ctg_data[5]
            result[ctg_from].append(
                __add_attributes(
                    [
                        f'{GFA_TO_KEY}:Z:{ctg_to}',
                        f'{GFA_FROM_ORIENT_KEY}:A:{ctg_from_orient}',
                        f'{GFA_TO_ORIENT_KEY}:A:{ctg_to_orient}',
                        f'{GFA_OVERLAP_KEY}:Z:{overlap}'
                    ],
                    [GFA_TO_KEY,GFA_FROM_ORIENT_KEY,GFA_TO_ORIENT_KEY,GFA_OVERLAP_KEY]
                )
            )
    return result

def write_GFA_to_FASTA(in_GFA_file, out_FASTA_file, in_gzipped, out_gzipped, sep=' '):
    """
    Create a FASTA file from a GFA file

    Args:
        - in_GFA_file (str): path to GFA file to read
        - out_FASTA_file (str): path to FASTA file to write
        - in_gzipped (bool): True if gzipped GFA file
        - out_gzipped (bool): True if FASTA file to be gzipped
        - sep (str): string for separating GFA attributes

    Returns:
        None, creates file out_FASTA_file
        header format: <contig name> <contig name>.GFA <attributes string=sep.join(key:value)>
        
    """    
    GFA_ctg_seqs = read_GFA_ctgs(
        in_GFA_file,
        attributes_list=['all'],
        gzipped=in_gzipped
    )
    ctg_records = [
        SeqRecord(
            Seq(y[GFA_SEQ_KEY]),
            id=x,
            name=x,
            description=f'{x}.GFA {__write_attributes(y, keys_to_remove=[GFA_SEQ_KEY])}'
        )
        for x,y in GFA_ctg_seqs.items()
    ]
    try:
        with __open_file_write(out_FASTA_file, gzipped=out_gzipped) as out_file:
            SeqIO.write(ctg_records, out_file, 'fasta')
    except Exception as e:
        process_exception(f'FASTA/GFA\tWriting {in_GFA_file} to {out_FASTA_file}: {e}')

