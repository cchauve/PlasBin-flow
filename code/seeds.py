import pandas as pd
from Bio import SeqIO
import numpy as np

def get_ctg_details(ctg):
    '''
    Parsing entry in the fasta file and getting contig details.
    Returns a dictionary with contig attributes as keys
    '''
    sequence = str(ctg.seq)
    gc_count = sequence.count('G') + sequence.count('C')
    length = len(sequence)
    return {
        'gc_count': gc_count,
        'length': length,
        'gc_percent': gc_count/length,
        'gd': 0,
        'type': 'chromosome'
    }

def parse_mapping(mapping_file):
    '''
    Computing the gene coverage intervals for each contig

    The mapfile is to be provided in BLAST output fmt 6. It is a tab separated file.
    Each row of the file has the following information in tab separated format
    qseqid      query or gene sequence id (str)
    sseqid      subject or contig sequence id (str)
    pident      percentage of identical positions (float)
    length      alignment or overlap length (int)
    mismatch    number of mismatches (int)
    gapopen     number of gap openings (int)
    qstart      start of alignment in query (int)
    qend        end of alignment in query (int)
    sstart      start of alignment in subject (int)
    send        end of alignment in subject (int)
    evalue      expect value (float)
    bitscore    bit score (int)
    '''
    covg_int = {}
    with open(mapping_file, 'r') as map:
        line = next(map)
        while line:
            tmp = line.split("\t")
            ctg = tmp[1]
            sstart, send = tmp[8], tmp[9]
            if ctg not in covg_int:
                covg_int[ctg] = []
            if int(sstart) > int(send):
                covg_int[ctg].append((int(send), int(sstart)))
            else:
                covg_int[ctg].append((int(sstart), int(send)))
            line = next(map, None)
    return covg_int

def get_union(intervals):
    '''
    Takes the gene covering intervals for a contig and finds their union
    The length of the union is used to compute gene coverage
    '''
    union = []
    for begin,end in sorted(intervals):
        if union and union[-1][1] >= begin-1:
            union[-1][1] = max(union[-1][1],end)
        else:
            union.append([begin,end])
    return union

def compute_gd(union, ctg_len):
    '''
    Computes gene density using list of coverage intervals and contig length
    '''
    covered = 0
    for interval in union:
        covered += interval[1] - interval[0] + 1
    return covered / ctg_len

#Classifying contigs as seeds according to gene density and length thresholds
def seed_by_param(gdt, lt, ctg_details):
    '''
    Function to classify contigs as seeds
    according to given threshold params (gene density and length)
    '''
    ctg_len, ctg_gd = ctg_details['length'], ctg_details['gd']
    return 1 if ctg_gd >= gdt and ctg_len >= lt else 0

def read_sample_assembly(all_ctgs_dict, sample, assembly_file):
    '''
    Reading assembly data
    Dictionary with key as the contig id
    Value as a dictionary of following attributes
    length (int), gd (float), gc_count (int), gc_percent (float), type (str)
    '''
    ctgs = SeqIO.parse(assembly_file,'fasta')
    for ctg in ctgs:
        ctg_id = f'{sample}_{ctg.id}'
        all_ctgs_dict[ctg_id] = get_ctg_details(ctg)
    return all_ctgs_dict

def get_gene_density(all_ctgs_dict, sample, gene_map_file):
    '''
    Reading gene to contig mapping information
    Computes gene density for all contigs
    '''
    covg_int = parse_mapping(gene_map_file)
    for ctg in covg_int:
        union = get_union(covg_int[ctg])
        ctg_id = f'{sample}_{ctg}'
        ctg_len = all_ctgs_dict[ctg_id]['length']
        all_ctgs_dict[ctg_id]['gd'] = compute_gd(union, ctg_len)
    return all_ctgs_dict

def get_ground_truth(all_ctgs_dict, all_pls_dict, sample, gt_file):
    '''
    Reading ground truth file
    Dictionary with key as the plamid id
    Value as a dictionary of following attributes
    sample (str), ctg_list (list)
    '''
    if gt_file != None:
        with open(gt_file, 'r') as gt:
            line = next(gt)
            while line:
                if line[0] != '#':
                    #Format assumption: (tab separated with first column PLS and second column CTG)
                    pls, ctg = line.split('\t')[0], line.split('\t')[1]
                    ctg_id = f'{sample}_{ctg}'
                    all_ctgs_dict[ctg_id]['type'] = 'plasmid'
                    try:
                        all_pls_dict[pls]['ctg_list'].append(ctg_id)
                    except KeyError:
                        all_pls_dict[pls] = {'ctg_list': [ctg_id], 'sample': sample}
                    line = next(gt, None)
    return all_ctgs_dict, all_pls_dict

#Reading reference sample files and storing data
def get_reference_data(input_csv_file):
    '''
    Takes csv file with addresses to assembly file,
    gene to contig mapping file (blast output) and ground truth file.
    Returns two dictionaries:
    1. Key: contig id, Value: nested dictionary with length, gd and contig source (plasmid/chromosome)
    2. Key: plasmid ids, Value: nested dictionary with sample name and list of contigs in the plasmid
    '''
    #Reading and storing input data for reference samples
    colnames = ['sample','assembly','mapping','ground_truth']
    PATHS_DF = pd.read_csv(input_csv_file)
    PATHS_DF.columns.values[[0, 1, 2, 3]] = colnames
    all_ctgs_dict = {}
    all_pls_dict = {}

    for index, row in PATHS_DF.iterrows():
        sample, assembly_file, gene_map_file, gt_file = row[0], row[1], row[2], row[3]
        all_ctgs_dict = read_sample_assembly(all_ctgs_dict, sample, assembly_file)
        all_ctgs_dict = get_gene_density(all_ctgs_dict, sample, gene_map_file)
        all_ctgs_dict, all_pls_dict = get_ground_truth(all_ctgs_dict, all_pls_dict, sample, gt_file)
    return all_ctgs_dict, all_pls_dict

def pls_seeds_by_thresholds(GD_THRESHOLDS, LEN_THRESHOLDS, all_ctgs_dict, all_pls_dict):
    '''
    Computing number of seeds for every plasmid
    Key: Plasmid ID
    Value: Nested dictionary with
    keys as theshold combinations (lt_gdt) and values as number of contigs classified as seeds using the thresholds
    '''
    seeds_dict = {}
    for pls in all_pls_dict:
        seeds_dict[pls] = {'sample': all_pls_dict[pls]['sample']}
        for gdt in GD_THRESHOLDS:
            for lt in LEN_THRESHOLDS:
                seeds_dict[pls][f'{lt}_{gdt}'] = 0
                for ctg_id in all_pls_dict[pls]['ctg_list']:
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
            for col_name, val in row.items():
                if val >= max_obj:
                    #best_params.add(row_name, col_name)
                    out.write(f'{row_name}\t{col_name}\n')

def compute_seeds_parameters_file(input_csv_file, out_file):
    #Reads reference data files and returns two dictionaries: one with contig details and another with plasmid details
    all_ctgs_dict, all_pls_dict = get_reference_data(input_csv_file)

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
