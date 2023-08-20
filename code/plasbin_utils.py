#!/usr/bin/env python

"""
Tuning/preprocessing for PlasBin-flow
Usage:

All commands take a common --log_file parameter (default=plasbin_utils.log).

Creating a database of plasmid genes from the samples plasmids
python plasbin_utils.py pls_genes_db --input_file input_file --out_dir out_dir --tmp_dir tmp_dir
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  pls_fasta: gzipped FASTA plasmids file
  where in FASTA files, the header of entries are GenBank accession
- out_dir: directory where the plasmid genes database written:
  pls.genes.fasta
- tmp_dir: temporary directory, not deleted

Mapping plasmid genes to samples contigs 
python plasbin_utils.py map_genes_to_ctgs --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --db_file pls_db_file
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- out_dir: directory where the mapping files are written
  <sample>.genes_mappings.txt
- tmp_dir: temporary directory, not deleted
- pls_db_file: path to plasmid genes database file

Computing ground truth for samples
python plasbin_utils.py ground_truth --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --out_file out_file --pid_threshold p --cov_threshold c
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  gfa: gzipped GFA file
  pls_fasta: gzipped FASTA plasmids file
- out_dir: directory where the ground truth files are written
  <sample>.ground_truth.tsv
- tmp_dir: temporary directory, not deleted
- out_file: path to new dataset file with ground truth files added
- p: percent identity threshold to define a mapping to a plasmid (default=0.95)
- c: coverage threshold to accept a blast hit  (default=0.8)

Computing GC content probabilities for samples
python plasbin_utils.py gc_probabilities --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --gc_intervals gc_intervals_file
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- out_dir: directory where the GC probabilities files are written
  <sample>.gc.tsv
- tmp_dir: temporary directory, not deleted
- gc_intervals_file: GC intervals file

Computing GC intervals from tuning samples
python plasbin_utils.py gc_intervals --input_file input_file --out_dir out_dir --tmp_dir tmp_dir
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  chr_fasta: gzipped chromosome FASTA file
  pls_fasta: gzipped FASTA plasmids file
- out_dir: directory where the GC intervals files are written
  gc.txt, gc.png
- tmp_dir: temporary directory, not deleted

Computing seeds parameters
python plasbin_utils.py seeds --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --db_file pls_db_file
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  gfa: gzipped GFA file
  ground_truth: ground truth file
- out_dir: directory where the seeds parameters file is written
  seeds.txt
- tmp_dir: temporary directory, not deleted
- pls_db_file: path to plasmid genes database file

python plasbin_tuning.py tuning --input_file input_file --out_dir out_dir --tmp_dir tmp_dir
- input_file: CSV file with one line per sample and 5 required fields:
  sample: sample name
  gfa: gzipped GFA file
  chr_fasta: gzipped chromosome FASTA file
  pls_fasta: gzipped FASTA plasmids file
  ground_truth: path to ground truth file
  where in FASTA files, the header of entries are GenBank accession
- output_dir: directory where the tuning files are written:
  pls.genes.fasta: plasmids genes database
  gc.txt: GC content per sample
  gc.png: GC content violin plot
  seeds.txt: seed parameters
- tmp_dir: temporary directory, not deleted
"""

"""
TODO: actually generate the GC intervals file
TODO: better logging and exceptions
TODO: clean temporary files
TODO: update README.md
"""

import sys
import os
import argparse
import shutil
import gzip
import subprocess
import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Reading input file

def read_samples(in_csv_file):
    """
    Reads the samples information CSV file
    
    Args:
        in_csv_file (str): path to the CSV file containing samples information

    Returns: 
        (DataFrame): dataframe indexed by sample id
    """
    try:
        samples_df = pd.read_csv(
            in_csv_file, sep = ',', header = 0, index_col = 'sample'
        )
        return samples_df
    except:
        logging.exception(
            f'Reading CSV dataset file {in_csv_file}'
        )
        sys.exit(1)

# Sample data access/modification functions
def _get_sample_col(samples_df, sample, col):
    try:
        return samples_df.at[sample,col]
    except:
        logging.exception(
            f'Reading {sample}.{col} in dataset file'
        )
        sys.exit(1)
def _get_gfa(samples_df, sample):
    return _get_sample_col(samples_df, sample, 'gfa')
def _get_chr_fasta(samples_df, sample):
    return _get_sample_col(samples_df, sample, 'chr_fasta')
def _get_pls_fasta(samples_df, sample):
    return _get_sample_col(samples_df, sample, 'pls_fasta')
def _get_ground_truth(samples_df, sample):
    return _get_sample_col(samples_df, sample, 'ground_truth')
def _set_ground_truth(samples_df, sample, gt_file):
    samples_df.at[sample,'ground_truth'] = gt_file

def _write_samples_df(samples_df, out_file):
    samples_df.to_csv(
        out_file, sep=',', header=True, index=True, index_label='sample'
    )

# Auxiliary functions

def _log_file(in_file):
    if os.path.isfile(in_file):
        logging.info(f'FILE\t{in_file}')
    else:
        logging.error(f'FILE\t{in_file} is missing')

def _convert_file(in_file, out_file, convert_fun, exception_msg):
    """
    Convert in_file into out_file, logging a warning if the file already exists

    Args:
       in_file (str): path to input file
       out_file (str): path to output file
       convert_fun (function): function of 2 parameters in_file, out_file
       exception_msg (str): message to log in case of an exception

    Returns:
       Creates out_file or log a warning if it already exists
    """
    if not os.path.isfile(out_file):
        try:
            convert_fun(in_file, out_file)
            _log_file(out_file)
        except:
            logging.exception(exception_msg)
            sys.exit(1)
    else:
        logging.warning(f'FILE {out_file} already exists')
        sys.exit(1)

def _gunzip_fasta(in_fasta_file, out_fasta_file):
    """
    Gunzip a FASTA file

    Args:
       in_fasta_file (str): path to input gzipped FASTA file
       out_fasta_file (str): path to output FASTA file

    Returns:
       Creates out_fasta_file or log a warning if it already exists
    """
    def _convert_fun(in_fasta_file, out_fasta_file):
        records = []
        with gzip.open(in_fasta_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                records.append(record)
            with open(out_fasta_file, 'w') as out_file:
                SeqIO.write(records, out_file, 'fasta')
    _convert_file(
        in_fasta_file, out_fasta_file, _convert_fun,
        f'Decompressing FASTA file {in_fasta_file} into {out_fasta_file}'
    )

def _gunzip_gfa(in_gfa_file, out_gfa_file):
    """
    Gunzip a GFA file

    Args:
       in_gfa_file (str): path to input gzipped GFA file
       out_gfa_file (str): path to output GFA file

    Returns:
       Creates out_gfa_file or log a warning if it already exists
    """
    def _convert_fun(in_gfa_file, out_gfa_file):
        with gzip.open(in_gfa_file) as in_file, open(out_gfa_file, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)
    _convert_file(
        in_gfa_file, out_gfa_file, _convert_fun,
        f'Decompressing GFA file {in_gfa_file} into {out_gfa_file}'
    )
        
def _gfa2fasta(in_gfa_file, out_fasta_file):
    """
    Convert an unzipped GFA file into an unzipped FASTA file

    Args:
       in_gfa_file (str): path to input gzipped GFA file
       out_fasta_file (str): path to output FASTA file

    Returns:
       Creates out_fasta_file or log a warning if it already exists
    """
    def _convert_fun(in_gfa_file, out_fasta_file):
        with open(in_gfa_file, 'r') as in_file, open(out_fasta_file, 'w') as out_file:
            for line in in_file.readlines():
                line_split = line.strip().split('\t')
                if line[0] == 'S':
                    out_file.write(
                        f'>{line_split[1]}\n{line_split[2]}\n'
                    )
    _convert_file(
        in_gfa_file, out_fasta_file, _convert_fun,
        f'Converting {in_gfa_file} to {out_fasta_file}'
    )

def _clean_files(files2clean):
    for in_file in files2clean:
        if os.path.isfile(in_file):
            os.remove(in_file)
        
def _create_directory(in_dir_list):
    for in_dir in in_dir_list:
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

def _run_cmd(cmd):
    cmd_str = ' '.join(cmd)
    logging.info(f'COMMAND {cmd_str}')
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info(f'STDOUT:\n{process.stdout}')
        logging.warning(f'STDERR:\n{process.stderr}')
    except subprocess.CalledProcessError as e:
        logging.exception(f'Running {cmd_str}')
        sys.exit(1)

def _read_seq_len(in_fasta_file):
    """
    Read contigs length from a FASTA file

    Args:
        in_fasta_file (str): path to (unzipped) FASTA file

    Returns:
        (Dictionary): contig id -> contig length
    """
    try:
        seq_len = {}
        with open(in_fasta_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                seq_len[record.id] = len(record.seq)
        return seq_len
    except:
        logging.exception(
            f'Reading FASTA file {in_fasta_file}'
        )
        sys.exit(1)

def _read_seq_id_gz(in_fasta_file):
    """
    Read contigs id from a gzipped FASTA file

    Args:
        in_fasta_file (str): path to gzipped FASTA file

    Returns:
        (List): contigs id
    """
    try:
        seq_id_list = []
        with gzip.open(in_fasta_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seq_id_list.append(record.id)
        return seq_id_list
    except:
        logging.exception(
            f'Reading gzipped FASTA file {in_fasta_file}'
        )
        sys.exit(1)

# File names

def _gfa_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.gfa')
def _gfa_fasta_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.gfa.fasta')
def _pls_fasta_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.pls.fasta')
def _pls_blastdb_prefix(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.pls.fasta.db')
def _pls_mappings_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.pls_mappings.txt')
def _pls_gb_file(in_dir):
    return os.path.join(in_dir, 'pls.genbank.txt')
def _genes_mappings_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.genes_mappings.txt')
def _ground_truth_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.ground_truth.tsv')
GC_FILE_PREFIX='gc'
def _gc_proba_file(in_dir, sample):
    return os.path.join(in_dir, f'{sample}.{GC_FILE_PREFIX}.tsv')
def _chr_pls_fasta_path_file(in_dir, file_type):
    return {
        'chr': os.path.join(in_dir, 'chr.fasta.txt'),
        'pls': os.path.join(in_dir, 'pls.fasta.txt')
    }[file_type]
def _seeds_input_file(in_dir):
    return os.path.join(in_dir, 'seeds_input.csv')
def _pls_genes_db_file(in_dir):
    return os.path.join(in_dir, 'pls.genes.fasta')
def _seeds_parameters_file(in_dir):
    return os.path.join(in_dir, 'seeds.txt')
def _gc_txt_file(in_dir):
    return os.path.join(in_dir, f'{GC_FILE_PREFIX}.txt')
def _gc_png_file(in_dir):
    return os.path.join(in_dir, f'{GC_FILE_PREFIX}.png')

# Processing functions

def create_tmp_data_files(tmp_dir, samples_df):
    """
    Creates a ungzipped GFA and FASTA file for each sample

    Args:
       samples_df (DataFrame): samples dataframe
       tmp_dir (str): path to temporary directory

    Returns:
       None, creates files _gfa_file(tmp_dir, sample) and _gfa_fasta_file(tmp_dir, sample) for each sample
    """
    logging.info(f'## Prepare temporary FASTA/GFA files')
    for sample in samples_df.index:
        logging.info(f'ACTION\tcopy assembly files for sample {sample}')
        gfa_file = _gfa_file(tmp_dir, sample)
        _gunzip_gfa(_get_gfa(samples_df,sample), gfa_file)
        _gfa2fasta(gfa_file, _gfa_fasta_file(tmp_dir, sample))

def _compute_ground_truth(out_dir, tmp_dir, sample, pid_threshold, cov_threshold):
    """
    Computes the ground truth file for a sample

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       sample (str): sample -> sample id
       pid_threshold (float): percent identity threshold
       cov_threshold (float): coverage threshold

    Returns:
      None, creates the file _ground_truth_file(out_dir, sample)
    """
    def _num_covered_positions(intervals):
        intervals.sort(key = lambda x: x[0])
        num_pos_covered,last_pos_covered = 0,0
        for start,end in intervals:
            if end <= last_pos_covered:
                pass
            else:
                num_pos_covered += end - max(last_pos_covered + 1, start) + 1
                last_pos_covered = end
        return num_pos_covered

    col_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore"
    ]  # outfmt 6
    hits = pd.read_csv(
        _pls_mappings_file(tmp_dir, sample),
        sep = '\t', names = col_names, dtype = str
    )
        
    ctg_len = _read_seq_len(_gfa_fasta_file(tmp_dir, sample))        
    pls_len = _read_seq_len(_pls_fasta_file(tmp_dir, sample))
    
    covered_sections = dict([(pred, []) for pred in ctg_len])	
    covered_per_ref = dict()		#Covered proportion per reference plasmid
    
    ref_ids = hits.sseqid.unique()	#Set of reference ids in the blast output
    
    for ref in sorted(ref_ids):
        covered_per_ref[ref] = {}
        for ctg in ctg_len.keys():
            covered_per_ref[ref][ctg] = []
            ctg_ref_hits = hits.loc[hits.qseqid == ctg].loc[hits.sseqid == ref]
            for index, row in ctg_ref_hits.iterrows():
                qstart, qend = int(row[6]), int(row[7])
                pident = float(row[2])/100
                interval = (qstart, qend) if qstart <= qend else (qend, qstart)
                if pident >= pid_threshold:
                    covered_sections[ctg].append(interval)
                    covered_per_ref[ref][ctg].append(interval)
    with open(_ground_truth_file(out_dir, sample), "w") as out_file:
        for ref in covered_per_ref.keys():
            for ctg in covered_per_ref[ref]:
                percent_mapping = _num_covered_positions(
                    covered_per_ref[ref][ctg]
                )/ctg_len[ctg]
                if percent_mapping >= cov_threshold:
                    pm = '{:.2f}'.format(percent_mapping)
                    out_file.write(
                        f'{ref}\t{ctg}\t{pm}\t{pls_len[ref]}\t{ctg_len[ctg]}\n'
                    )

def create_ground_truth_files(
        out_dir, tmp_dir, samples_df, out_file,
        pid_threshold=0.95, cov_threshold=0.8
):
    """
    Creates a TSV ground truth file per sample

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       out_file (str): name of dataset file augmeted with ground truth files
       pid_threshold (float): percent identity threshold to define a true psitive mapping
       cov_threshold (float): coverage threshold to consider a hit 

    Returns:
       None, creates files _ground_truth_file(out_dir, sample) for each sample  
    """
    logging.info(f'## Compute ground truth and create new dataset CSV file')
    for sample in samples_df.index:
        logging.info(f'ACTION\tground truth for {sample}')
        pls_fasta_file = _pls_fasta_file(tmp_dir, sample)
        _gunzip_fasta(_get_pls_fasta(samples_df, sample), pls_fasta_file)
        logging.info(f'ACTION\tcompute blast database for {pls_fasta_file}')
        pls_blastdb_prefix = _pls_blastdb_prefix(tmp_dir, sample)
        gfa_fasta_file = _gfa_fasta_file(tmp_dir, sample)
        pls_mappings_file = _pls_mappings_file(tmp_dir, sample)
        cmd1 = [
            'makeblastdb',
            '-in', pls_fasta_file,
            '-dbtype', 'nucl',
            '-out', pls_blastdb_prefix
        ]
        _run_cmd(cmd1)
        logging.info(f'ACTION\tmap {gfa_fasta_file} to {pls_blastdb_prefix}')        
        cmd2 = [
            'blastn', '-task', 'megablast',
            '-query', gfa_fasta_file,
            '-db', pls_blastdb_prefix,
            '-out', pls_mappings_file,
            '-outfmt', '6'
        ]
        _run_cmd(cmd2)
        _log_file(pls_mappings_file)
        logging.info(f'ACTION\tcompute ground truth file')                
        _compute_ground_truth(
            out_dir, tmp_dir, sample, pid_threshold, cov_threshold
        )
        ground_truth_file = _ground_truth_file(out_dir, sample)
        _set_ground_truth(samples_df, sample, ground_truth_file)
        _log_file(ground_truth_file)
    _write_samples_df(samples_df, out_file)
    _log_file(out_file)
    
def create_pls_genes_db(out_dir, tmp_dir, samples_df):
    """
    Creates a plasmid genes database

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe

    Returns:
       None, creates file _pls_genes_db_file(out_dir, out_file_name)
    """
    def _create_input_file(samples_df, input_file):
        try:
            with open(input_file, 'w') as out_file:
                for sample in samples_df.index:
                    for seq_id in _read_seq_id_gz(
                            _get_pls_fasta(samples_df, sample)
                    ):
                        out_file.write(f'{seq_id}\n')
            _log_file(input_file)
        except:
            logging.exception(f'Creating {input_file}')
            sys.exit(1)
    
    logging.info(f'## Compute plasmid genes database')
    logging.info(f'ACTION\tcreate plasmid GenBank accessions file')
    pls_gb_file = _pls_gb_file(tmp_dir)
    _create_input_file(samples_df, pls_gb_file)
    logging.info(f'ACTION\tprocess {pls_gb_file}')
    pls_genes_db_file = _pls_genes_db_file(out_dir)
    cmd = [
        'python', 'get_gd.py', 'create',
        pls_genes_db_file,
        '-a', pls_gb_file
    ]
    _run_cmd(cmd)
    _log_file(pls_genes_db_file)

def map_pls_genes_to_contigs(out_dir, tmp_dir, samples_df, db_file):
    """
    Map plasmid genes to samples contigs

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       db_file (str): path to plasmid genes database file

    Returns:
       None, creates files _genes_mappings_file(out_dir, sample)nfor each sample
    """
    logging.info(f'## Mapping sample contigs to plasmid genes database')
    for sample in samples_df.index:
        genes_mappings_file = _genes_mappings_file(out_dir, sample)
        logging.info(f'ACTION\tmapping {sample} to {db_file}')
        cmd = [
            'python', 'get_gd.py', 'map',
            db_file,
            genes_mappings_file,
            '-a', _gfa_file(tmp_dir, sample)
        ]
        _run_cmd(cmd)
        _log_file(genes_mappings_file)
        
def create_GC_content_intervals_file(out_dir, tmp_dir, samples_df):
    """
    Creates GC content intervals files

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       out_file_prefix (str): prefix of output file names

    Returns:
       None, creates files _gc_txt_file(out_dir) and _gc_png_file(out_dir)
    """

    def _create_input_file(samples_df, input_file, file_type):
        logging.info(f'ACTION\trecord {file_type} FASTA files paths')
        try:
            __get_file = {'chr': _get_chr_fasta, 'pls': _get_pls_fasta}
            with open(input_file, 'w') as out_file:
                for sample in samples_df.index:
                    out_file.write(
                        f'{__get_file[file_type](samples_df, sample)}\n'
                    )
            _log_file(input_file)
        except:
            logging.exception(f'Creating {input_file}')
            sys.exit(1)
    
    logging.info(f'## Compute GC content intervals files')
    for file_type in ['chr','pls']:
        input_file = _chr_pls_fasta_path_file(tmp_dir, file_type)
        _create_input_file(samples_df, input_file, file_type)
    logging.info(f'ACTION\tcompute GC content intervals files')
    out_txt_file = _gc_txt_file(out_dir)
    out_png_file = _gc_png_file(out_dir)
    cmd = [
        'python', 'analyse_GC_content.py',
        '--chr', _chr_pls_fasta_path_file(tmp_dir, 'chr'),
        '--pls', _chr_pls_fasta_path_file(tmp_dir, 'pls'),        
        '--out', out_txt_file,
        '--vplot', out_png_file
    ]
    _run_cmd(cmd)
    _log_file(out_txt_file)
    _log_file(out_png_file)

def create_GC_content_probabilities_file(
    out_dir, tmp_dir, gc_intervals_file, samples_df
): 
    """
    Creates GC content probabilities files

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe

    Returns:
       None, creates files _gc_proba_file(out_dir, sample) for each sample
    """
    logging.info(f'## GC content probabilities files')
    for sample in samples_df.index:
        logging.info(f'ACTION\tcompute GC content probabilities file for sample {sample}')
        gc_proba_file = _gc_proba_file(out_dir, sample)
        cmd = [
            'python', 'get_gc_probs.py',
            '-ag', _gfa_file(tmp_dir, sample),
            '-outfile', gc_proba_file,
            '-gcint', gc_intervals_file
        ]
        _run_cmd(cmd)
        _log_file(gc_proba_file)

def create_seeds_parameters_file(out_dir, tmp_dir, samples_df, db_file):
    """
    Creates a file containing the parameters defining seeds

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       db_file (str): path to plasmid genes database file

    Returns:
       None, creates file _seeds_parameters_file(out_dir, out_file_name)
    """
    def _create_input_file(samples_df, input_file):
        try:
            with open(input_file, 'w') as out_file:
                for sample in samples_df.index:
                    mappings_file = _genes_mappings_file(tmp_dir, sample)
                    gt_file = _get_ground_truth(samples_df, sample)
                    fasta_file = _gfa_fasta_file(tmp_dir, sample)
                    out_file.write(
                        f'{sample},{fasta_file},{mappings_file},{gt_file}\n'
                    )
            _log_file(input_file)
        except:
            logging.exception(f'Creating {input_file}')
            sys.exit(1)
    
    logging.info(f'## Compute seeds parameters file')
    logging.info(f'ACTION\tmap plasmid genes to contigs')    
    map_pls_genes_to_contigs(
        tmp_dir, tmp_dir, samples_df, db_file
    )
    logging.info(f'ACTION\tcreate plasmids seeds input file')    
    seeds_input_file = _seeds_input_file(tmp_dir)
    _create_input_file(samples_df, seeds_input_file)
    logging.info(f'ACTION\tcreate seeds parameters file') 
    seeds_parameters_file = _seeds_parameters_file(out_dir)
    cmd = [
        'python', 'analyse_seed_eligibility.py',
        '--paths', seeds_input_file,
        '--out', seeds_parameters_file
    ]
    _run_cmd(cmd)
    _log_file(seeds_parameters_file)


def main():
    argparser = argparse.ArgumentParser(description='PlasBin-flow utils')
    argparser.add_argument('--input_file', type=str, help='Samples CSV file')
    argparser.add_argument('--out_dir', type=str, help='Output directory')    
    argparser.add_argument('--tmp_dir', type=str, help='Temporary directory')
    argparser.add_argument('--log_file', type=str, default='plasbin_utils.log', help='Log file')
    subparsers = argparser.add_subparsers(title='commands')
    # Creating a plasmid genes database
    db_parser = subparsers.add_parser('pls_genes_db', parents=[argparser], add_help=False)
    db_parser.set_defaults(cmd='pls_genes_db')
    # Map plasmid genes to samples contigs
    g2c_parser = subparsers.add_parser('map_genes_to_ctgs', parents=[argparser], add_help=False)
    g2c_parser.set_defaults(cmd='map_genes_to_ctgs')
    g2c_parser.add_argument('--db_file', type=str, help='Plasmids genes database FASTA file')
    # Computing ground truth files
    gt_parser = subparsers.add_parser('ground_truth', parents=[argparser], add_help=False)
    gt_parser.set_defaults(cmd='ground_truth')
    gt_parser.add_argument('--out_file', type=str, help='Path to dataset file with added ground truth files')    
    gt_parser.add_argument('--pid_threshold', type=float, default=0.95, help='Percent identity threshold in [0,1]')
    gt_parser.add_argument('--cov_threshold', type=float, default=0.8, help='Percent coverage threshold in [0,1]')    
    # Computing seeds parameters
    seeds_parser = subparsers.add_parser('seeds', parents=[argparser], add_help=False)
    seeds_parser.set_defaults(cmd='seeds')
    seeds_parser.add_argument('--db_file', type=str, help='Plasmids genes database FASTA file')
    # Computing GC contents intervals
    gci_parser = subparsers.add_parser('gc_intervals', parents=[argparser], add_help=False)
    gci_parser.set_defaults(cmd='gc_intervals')
    # Computing GC contents probabilities
    gcp_parser = subparsers.add_parser('gc_probabilities', parents=[argparser], add_help=False)
    gcp_parser.set_defaults(cmd='gc_probabilities')
    gcp_parser.add_argument('--gc_intervals', type=str, help='GC content intervals file')
    # Tuning
    tuning_parser = subparsers.add_parser('tuning', parents=[argparser], add_help=False)
    tuning_parser.set_defaults(cmd='tuning')

    args = argparser.parse_args()

    logging.basicConfig(
        filename=args.log_file,
        filemode='w',
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s'
    )

    if args.cmd == 'pls_genes_db':
        samples_df = read_samples(args.input_file)
        files2clean = [
            _pls_genes_db_file(args.out_dir)
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_pls_genes_db(
            args.out_dir, args.tmp_dir, samples_df
        )

    elif args.cmd == 'map_genes_to_ctgs':
        samples_df = read_samples(args.input_file)
        files2clean = [
            _genes_mappings_file(args.out_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        map_pls_genes_to_contigs(
            args.out_dir, args.tmp_dir, samples_df, args.db_file
        )

    elif args.cmd == 'ground_truth':
        samples_df = read_samples(args.input_file)
        files2clean = [
            _ground_truth_file(args.out_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        create_ground_truth_files(
            args.out_dir, args.tmp_dir, samples_df, args.out_file
        )

    elif args.cmd == 'seeds':
        samples_df = read_samples(args.input_file)
        files2clean = [_seeds_parameters_file(args.out_dir)]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        create_seeds_parameters_file(
            args.out_dir, args.tmp_dir, samples_df,
            args.db_file
        )

    elif args.cmd == 'gc_intervals':
        samples_df = read_samples(args.input_file)
        files2clean = [
            _gc_txt_file(args.out_dir),
            _gc_png_file(args.out_dir)
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_GC_content_intervals_file(
            args.out_dir, args.tmp_dir, samples_df
        )

    elif args.cmd == 'gc_probabilities':
        samples_df = read_samples(args.input_file)
        files2clean = [
            _gc_proba_file(args.out_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_GC_content_probabilities_file(
            args.out_dir, args.tmp_dir,
            args.gc_intervals, samples_df
        )
        
    elif args.cmd == 'tuning':
        samples_df = read_samples(args.input_file)
        files2clean = [
            _pls_genes_db_file(args.out_dir),
            _seeds_parameters_file(args.out_dir),
            _gc_txt_file(args.out_dir),
            _gc_png_file(args.out_dir)
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        create_pls_genes_db(
            args.out_dir, args.tmp_dir, samples_df
        )
        create_seeds_parameters_file(
            args.out_dir, args.tmp_dir, samples_df,
            _pls_genes_db_file(args.out_dir)
        )
        create_GC_content_intervals_file(
            args.out_dir, args.tmp_dir, samples_df
        )
        
if __name__ == "__main__":
    main()
