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
- out_file: [optional] path to new dataset file with mappings files added

Computing ground truth for samples
python plasbin_utils.py ground_truth --input_file input_file --out_dir out_dir --tmp_dir tmp_dir [--out_file out_file --pid_threshold p --cov_threshold c]
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  gfa: gzipped GFA file
  pls_fasta: gzipped FASTA plasmids file
- out_dir: directory where the ground truth files are written
  <sample>.ground_truth.tsv
- tmp_dir: temporary directory, not deleted
- out_file: [optional] path to new dataset file with ground truth files added
- p: [optional] percent identity threshold to define a mapping to a plasmid (default=0.95)
- c: [optional] coverage threshold to accept a blast hit  (default=0.8)

Computing GC content probabilities for samples
python plasbin_utils.py gc_probabilities --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --gc_intervals gc_intervals_file
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- out_dir: directory where the GC probabilities files are written
  <sample>.gc.tsv
- tmp_dir: temporary directory, not deleted
- out_file: [optional] path to new dataset file with ground truth files added
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
- n_gcints: number of GC content intervals between 0 and 1, default value: 6

Computing gene density
python plasbin_utils.py gene_density --input_file input_file --out_dir out_dir --tmp_dir tmp_dir [--out_file out_file]
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  gfa: gzipped GFA file
  gene_mappings: plasmid genes to contigs mappings file
- out_dir: directory where the gene density files are written
  (one per sample, <sample>.gene_density.tsv)
- tmp_dir: temporary directory, not deleted
- out_file: [optional] path to new dataset file with gene density files added

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
- n_gcints: number of GC content intervals between 0 and 1, default value: 6

python plasbin_tuning.py preprocessing --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --pls_db pls_db --gc_intervals gc_intervals --out_file out_file
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- output_dir: directory where the tuning files are written:
  plasmid genes to contigs mappings (one per sample, <sample>.genes_mappings.txt)
  GC content probabilities files  (one per sample, <sample>.gc.tsv)
- tmp_dir: temporary directory, not deleted
- pls_db: plasmid genes database
- gc_intervals: GC content intervals file
- out_file: augmented dataset CSV file, with mappings and GC probabilities files added
"""

import sys
import os
import argparse
import shutil
import gzip
import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import create_db as cd

from gc_content import (
    compute_gc_intervals_files,
    compute_gc_probabilities_file
)
from ground_truth import (
    compute_ground_truth_file
)
from gene_density import (
    compute_gene_density_file
)
from seeds import (
    compute_seeds_parameters_file
)
from gfa_fasta_utils import (
    gunzip_GFA,
    gunzip_FASTA,
    write_GFA_to_FASTA,
    read_FASTA_len,
    read_FASTA_id
)
from mappings_utils import run_blast6
from log_errors_utils import (
    check_file,
    log_file,
    CustomException,
    clean_files,
    create_directory,
    run_cmd,
    process_error,
    process_exception
)

# Reading input file

GFA_COL = 'gfa'
CHR_COL = 'chr_fasta'
PLS_COL = 'pls_fasta'
GT_COL = 'ground_truth'
GC_COL = 'gc_probabilities'
MAPPINGS_COL = 'genes2ctgs_mappings'
GD_COL = 'gene_density'

def check_input_files(samples_df):
    """
    Check that all file entries in sample_df do exist
    Args:
        - sample_df (DataFrame)
    Returns:
        - List(str): empty entries
        - List(str): missing files
    """
    empty_files = []
    missing_files = []
    for sample,files in samples_df.iterrows():
        for col,file_col in files.items():
            if pd.isnull(file_col):
                empty_files.append(f'{sample}.{col}')
            elif not os.path.isfile(file_col):
                missing_files.append(file_col)
    return empty_files,missing_files

def read_samples(in_csv_file, required_columns):
    """
    Reads the samples information CSV file
    
    Args:
        in_csv_file (str): path to the CSV file containing samples information

    Returns: 
        (DataFrame): dataframe indexed by sample id
    """
    try:
        samples_df = pd.read_csv(
            in_csv_file,
            sep = ',',
            header = 0,
            index_col = 'sample',
        )
        missing_columns = [
            col for col in required_columns
            if col not in samples_df.columns
        ]
        if len(missing_columns) > 0:
            msg = ' '.join(missing_columns)
            raise CustomException(f'{in_csv_file}: Missing column(s) {msg}') 
        empty_files,missing_files = check_input_files(samples_df)
        if len(empty_files) > 0:
            msg = ' '.join(empty_files)
            raise CustomException(f'{in_csv_file}: Empty file entry {msg}')
        if len(missing_files) > 0:
            msg = ' '.join(missing_files)
            raise CustomException(f'{in_csv_file}: Missing files {msg} ')
    except Exception as e:
        process_exception(f'Reading CSV dataset file {in_csv_file}: {e}')
    else:
        return samples_df
        
# Sample data access/modification functions
def _get_sample_col(samples_df, sample, col):
    """
    Args:
        samples_df (DataFrame): samples dataframe
        sample (str): name of a sample
        col (str): name of a column of the dataframe

    Returns:
        Value of field col for sample
    """
    return samples_df.at[sample,col]
        
def _get_gfa(samples_df, sample):
    """ Path to gzipped GFA file """
    return _get_sample_col(samples_df, sample, GFA_COL)
def _get_chr_fasta(samples_df, sample):
    """ Path to gzipped chromosome FASTA file """
    return _get_sample_col(samples_df, sample, CHR_COL)
def _get_pls_fasta(samples_df, sample):
    """ Path to gzipped plasmids FASTA file """
    return _get_sample_col(samples_df, sample, PLS_COL)
def _get_ground_truth(samples_df, sample):
    """ Path to ground truth file """
    return _get_sample_col(samples_df, sample, GT_COL)
def _set_ground_truth(samples_df, sample, gt_file):
    """ Set path to ground truth file for sample """
    samples_df.at[sample,GT_COL] = gt_file
def _get_gc_prob(samples_df, sample):
    """ Path to GC content probabilities file """
    return _get_sample_col(samples_df, sample, GC_COL)
def _set_gc_prob(samples_df, sample, gc_prob_file):
    """ Set path to GC content probabilities file for sample """
    samples_df.at[sample,GC_COL] = gc_prob_file
def _get_genes2ctgs_mappings(samples_df, sample):
    """ Path to genes to contigs mappings file """
    return _get_sample_col(samples_df, sample, MAPPINGS_COL)
def _set_genes2ctgs_prob(samples_df, sample, mappings_file):
    """ Set path to genes to contigs mappings file """
    samples_df.at[sample,MAPPINGS_COL] = mappings_file
def _get_gene_density(samples_df, sample):
    """ Path to gene density file """
    return _get_sample_col(samples_df, sample, GD_COL)
def _set_gene_density(samples_df, sample, gd_file):
    """ Set path to gene density file """
    samples_df.at[sample,GD_COL] = gd_file

    
def _write_samples_df(samples_df, out_file):
    """
    Write samples DataFrame samples_df to CSV file out_file
    """
    samples_df.to_csv(
        out_file, sep=',', header=True, index=True, index_label='sample'
    )

# File names

## Samples specific file path
def _gfa_file(in_dir, sample):
    """ Unzipped GFA file """
    return os.path.join(in_dir, f'{sample}.gfa')
def _gfa_fasta_file(in_dir, sample):
    """ Unzipped FASTA file for contigs in the GFA file """
    return os.path.join(in_dir, f'{sample}.gfa.fasta')
def _pls_fasta_file(in_dir, sample):
    """ Gzipped plasmids FASTA file """
    return os.path.join(in_dir, f'{sample}.pls.fasta')
def _pls_mappings_file(in_dir, sample):
    """ Mapping of sample contigs to plasmids """
    return os.path.join(in_dir, f'{sample}.pls_mappings.txt')
def _genes_mappings_file(in_dir, sample):
    """ Mapping of genes in plasmid genes database to sample contigs """
    return os.path.join(in_dir, f'{sample}.genes_mappings.txt')
def _ground_truth_file(in_dir, sample):
    """ Ground truth file """
    return os.path.join(in_dir, f'{sample}.ground_truth.tsv')
def _gene_density_file(in_dir, sample):
    """ Gene density file """
    return os.path.join(in_dir, f'{sample}.gd.tsv')
GC_FILE_PREFIX='gc'
def _gc_proba_file(in_dir, sample):
    """ GC probabilities file """
    return os.path.join(in_dir, f'{sample}.{GC_FILE_PREFIX}.tsv')
## Datasets wide files
def _pls_genes_db_file(in_dir):
    """ Plasmids genes FASTA database file """
    return os.path.join(in_dir, 'pls.genes.fasta')
def _seeds_parameters_file(in_dir):
    """ Seeds parameters file """
    return os.path.join(in_dir, 'seeds.txt')
def _gc_csv_file(in_dir):
    """ GC content intervals TXT file """
    return os.path.join(in_dir, f'{GC_FILE_PREFIX}.csv')
def _gc_png_file(in_dir):
    """ GC content violin plot PNG file """
    return os.path.join(in_dir, f'{GC_FILE_PREFIX}.png')
def _gc_intervals_file(in_dir):
    """ GC content intervals TXT file """
    return os.path.join(in_dir, f'{GC_FILE_PREFIX}.txt')

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
        gunzip_GFA(_get_gfa(samples_df,sample), gfa_file)
        log_file(gfa_file)
        gfa_fasta_file = _gfa_fasta_file(tmp_dir, sample)
        write_GFA_to_FASTA(
            gfa_file, gfa_fasta_file,
            in_gzipped=False, out_gzipped=False
        )
        log_file(gfa_fasta_file)

def create_ground_truth_files(
        out_dir, tmp_dir, samples_df,
        pid_threshold=0.95, cov_threshold=0.8
):
    """
    Creates a TSV ground truth file per sample

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       pid_threshold (float): percent identity threshold to define a true positive mapping
       cov_threshold (float): coverage threshold to consider a hit 

    Returns:
       None, creates files _ground_truth_file(out_dir, sample) for each sample  

    Assumptions:
       makeblastdb and blastn are in the default path
    """
    logging.info(f'## Compute ground truth and create new dataset CSV file')
    for sample in samples_df.index:
        logging.info(f'ACTION\tground truth for {sample}')
        logging.info(f'ACTION\tMapping contigs to plasmids')
        pls_fasta_file = _pls_fasta_file(tmp_dir, sample)
        gunzip_FASTA(_get_pls_fasta(samples_df, sample), pls_fasta_file)
        gfa_fasta_file = _gfa_fasta_file(tmp_dir, sample)
        pls_mappings_file = _pls_mappings_file(tmp_dir, sample)
        run_blast6(gfa_fasta_file, pls_fasta_file, pls_mappings_file)
        logging.info(f'ACTION\tcompute ground truth file')                
        ground_truth_file = _ground_truth_file(out_dir, sample)
        compute_ground_truth_file(
            sample,
            pls_mappings_file,
            read_FASTA_len(_gfa_fasta_file(tmp_dir, sample), gzipped=False),
            read_FASTA_len(_pls_fasta_file(tmp_dir, sample), gzipped=False),
            pid_threshold, cov_threshold,
            ground_truth_file 
        )
        _set_ground_truth(samples_df, sample, ground_truth_file)
        log_file(ground_truth_file)

def create_gene_density_files(out_dir, samples_df):
    """
    Creates a TSV gene density file per sample

    Args:
       out_dir (str): path to output directory
       samples_df (DataFrame): samples dataframe

    Returns:
       None, creates files _gene_density_file(out_dir, sample) for each sample  

    """
    logging.info(f'## Compute gene densit and create new dataset CSV file')
    for sample in samples_df.index:
        logging.info(f'ACTION\tcompute gene density for {sample}')
        gfa_file = _get_gfa(samples_df, sample)
        mappings_file = _get_genes2ctgs_mappings(samples_df, sample)
        gd_out_file = _gene_density_file(out_dir, sample)
        compute_gene_density_file(
            gfa_file, mappings_file, gd_out_file, gfa_gzipped=True
        )
        _set_gene_density(samples_df, sample, gd_out_file)
        log_file(gd_out_file)
    
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
        with open(input_file, 'w') as out_file:
            for sample in samples_df.index:
                pls_fasta_file = _get_pls_fasta(samples_df, sample)
                try:
                    pls_ids = read_FASTA_id(pls_fasta_file, gzipped=True)
                except Exception as e:
                    process_exception(
                        f'Reading plasmids FASTA file {pls_fasta_file}: {e}'
                    )                    
                else:
                    pls_ids_str = '\n'.join(pls_ids)
                    out_file.write(f'{pls_ids_str}\n')
        log_file(input_file)

    logging.info(f'## Compute plasmid genes database')
    logging.info(f'ACTION\tcreate plasmid GenBank accessions file')
    pls_gb_file = os.path.join(tmp_dir, 'pls.genbank.txt')
    _create_input_file(samples_df, pls_gb_file)
    logging.info(f'ACTION\tprocess {pls_gb_file}')
    pls_genes_db_file = _pls_genes_db_file(out_dir)
    cd.create(
        pls_genes_db_file,
        pls_gb_file)
    log_file(pls_genes_db_file)

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
    logging.info(f'## Mapping plasmid genes database to sample contigs')
    for sample in samples_df.index:
        logging.info(f'ACTION Computing genes to contigs mappings for {sample}')
        genes_mappings_file = _genes_mappings_file(out_dir, sample)
        gfa_fasta_file = _gfa_fasta_file(tmp_dir, sample)
        run_blast6(db_file, gfa_fasta_file, genes_mappings_file)
        _set_genes2ctgs_prob(samples_df, sample, genes_mappings_file)
        
def create_GC_content_intervals_file(out_dir, tmp_dir, samples_df, n_gcints):
    """
    Creates GC content intervals files

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       out_file_prefix (str): prefix of output file names
       n_gcints (int): number of intervals

    Returns:
       None, creates files _gc_csv_file(out_dir) and _gc_png_file(out_dir)
    """

    def _chr_pls_fasta_path_file(in_dir, file_type):
        """ List of paths to chromosomes/plasmids gzipped FASTA files """
        return {
            'chr': os.path.join(in_dir, 'chr.fasta.txt'),
            'pls': os.path.join(in_dir, 'pls.fasta.txt')
        }[file_type]
    
    def _create_input_file(samples_df, input_file, file_type):
        logging.info(f'ACTION\trecord {file_type} FASTA files paths')
        __get_file = {'chr': _get_chr_fasta, 'pls': _get_pls_fasta}
        with open(input_file, 'w') as out_file:
            for sample in samples_df.index:
                out_file.write(
                    f'{__get_file[file_type](samples_df, sample)}\n'
                )
        log_file(input_file)
    
    logging.info(f'## Compute GC content intervals files')
    for file_type in ['chr','pls']:
        input_file = _chr_pls_fasta_path_file(tmp_dir, file_type)
        _create_input_file(samples_df, input_file, file_type)
    out_csv_file = _gc_csv_file(out_dir)
    out_png_file = _gc_png_file(out_dir)
    out_intervals_file = _gc_intervals_file(out_dir)
    logging.info(f'ACTION\tcompute GC content intervals files {out_intervals_file}')
    compute_gc_intervals_files(
        _chr_pls_fasta_path_file(tmp_dir, 'chr'),
        _chr_pls_fasta_path_file(tmp_dir, 'pls'),
        out_csv_file,
        out_png_file,
        out_intervals_file,
        n_gcints
    )
    log_file(out_csv_file)
    log_file(out_png_file)
    log_file(out_intervals_file)

def create_GC_content_probabilities_files(
        out_dir, tmp_dir, gc_intervals_file, samples_df
): 
    """
    Creates GC content probabilities files

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       gc_intervals_file (str): path to file with GC content intervals
       samples_df (DataFrame): samples dataframe

    Returns:
       None, creates files _gc_proba_file(out_dir, sample) for each sample
    """
    logging.info(f'## GC content probabilities files')
    for sample in samples_df.index:
        logging.info(f'ACTION\tcompute GC content probabilities file for sample {sample}')
        gc_proba_file = _gc_proba_file(out_dir, sample)
        compute_gc_probabilities_file(
            _gfa_file(tmp_dir, sample),
            gc_intervals_file,
            gc_proba_file            
        )
        _set_gc_prob(samples_df, sample, gc_proba_file)        
        log_file(gc_proba_file)
    
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
        with open(input_file, 'w') as out_file:
            for sample in samples_df.index:
                mappings_file = _genes_mappings_file(tmp_dir, sample)
                gt_file = _get_ground_truth(samples_df, sample)
                fasta_file = _gfa_fasta_file(tmp_dir, sample)
                out_file.write(
                    f'{sample},{fasta_file},{mappings_file},{gt_file}\n'
                )
        log_file(input_file)
    
    logging.info(f'## Compute seeds parameters file')
    logging.info(f'ACTION\tmap plasmid genes to contigs')    
    map_pls_genes_to_contigs(
        tmp_dir, tmp_dir, samples_df, db_file
    )
    logging.info(f'ACTION\tcreate plasmids seeds input file')    
    seeds_input_file = os.path.join(tmp_dir, 'seeds_input.csv')
    _create_input_file(samples_df, seeds_input_file)
    logging.info(f'ACTION\tcreate seeds parameters file') 
    seeds_parameters_file = _seeds_parameters_file(out_dir)
    compute_seeds_parameters_file(seeds_input_file, seeds_parameters_file)
    log_file(seeds_parameters_file)


# Main

def _read_input(cmd, input_file):
    check_file(input_file)
    required_columns = {
        'pls_genes_db': [PLS_COL],
        'map_genes_to_ctgs': [GFA_COL],
        'ground_truth': [GFA_COL,PLS_COL],
        'gene_density': [GFA_COL,MAPPINGS_COL],
        'seeds': [GFA_COL,GT_COL],
        'gc_intervals': [PLS_COL,CHR_COL],
        'gc_probabilities': [GFA_COL],
        'tuning': [GFA_COL,PLS_COL,CHR_COL,GT_COL],
        'preprocessing': [GFA_COL]
    }
    samples_df = read_samples(input_file, required_columns[cmd])
    return samples_df

def _clean_output_files(cmd, samples_df, out_dir, tmp_dir):
    samples_list = samples_df.index
    files2clean = {
        'pls_genes_db': [_pls_genes_db_file(out_dir)],
        'map_genes_to_ctgs': [
            _genes_mappings_file(out_dir, sample)
            for sample in samples_list
        ],
        'ground_truth': [
            _ground_truth_file(out_dir, sample)
            for sample in samples_list
        ],
        'gene_density': [
            _gene_density_file(out_dir, sample)
            for sample in samples_list
        ],
        'seeds': [_seeds_parameters_file(out_dir)],
        'gc_intervals': [
            _gc_csv_file(out_dir),
            _gc_png_file(out_dir),
            _gc_intervals_file(out_dir)
        ],
        'gc_probabilities': [
            _gc_proba_file(out_dir, sample)
            for sample in samples_list
        ],
        'tuning': [
            _pls_genes_db_file(out_dir),
            _seeds_parameters_file(out_dir),
            _gc_csv_file(out_dir),
            _gc_png_file(out_dir),
            _gc_intervals_file(out_dir)
        ],
        'preprocessing': [
            _genes_mappings_file(out_dir, sample)
            for sample in samples_list
        ] + [
            _gc_proba_file(out_dir, sample)
            for sample in samples_list
        ] + [
             _gene_density_file(out_dir, sample)
            for sample in samples_list
        ]
    }
    clean_files(files2clean[cmd])

def _read_arguments():
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
    g2c_parser.add_argument('--out_file', type=str, help='Path to dataset file with added mappings files')        
    # Computing ground truth files
    gt_parser = subparsers.add_parser('ground_truth', parents=[argparser], add_help=False)
    gt_parser.set_defaults(cmd='ground_truth')
    gt_parser.add_argument('--out_file', type=str, help='Path to dataset file with added ground truth files')    
    gt_parser.add_argument('--pid_threshold', type=float, default=0.95, help='Percent identity threshold in [0,1]')
    gt_parser.add_argument('--cov_threshold', type=float, default=0.8, help='Percent coverage threshold in [0,1]')    
    # Computing gene density files
    gd_parser = subparsers.add_parser('gene_density', parents=[argparser], add_help=False)
    gd_parser.set_defaults(cmd='gene_density')
    gd_parser.add_argument('--out_file', type=str, help='Path to dataset file with added gene density files')    
    # Computing seeds parameters
    seeds_parser = subparsers.add_parser('seeds', parents=[argparser], add_help=False)
    seeds_parser.set_defaults(cmd='seeds')
    seeds_parser.add_argument('--db_file', type=str, help='Plasmids genes database FASTA file')
    # Computing GC contents intervals
    gci_parser = subparsers.add_parser('gc_intervals', parents=[argparser], add_help=False)
    gci_parser.set_defaults(cmd='gc_intervals')
    gci_parser.add_argument('--n_gcints', type=int, default=6, help='Number of GC content intervals between 0 and 1')
    # Computing GC contents probabilities
    gcp_parser = subparsers.add_parser('gc_probabilities', parents=[argparser], add_help=False)
    gcp_parser.set_defaults(cmd='gc_probabilities')
    gcp_parser.add_argument('--out_file', type=str, help='Path to dataset file with added GC probabilities files')        
    gcp_parser.add_argument('--gc_intervals', type=str, help='GC content intervals file')
    # Tuning
    tuning_parser = subparsers.add_parser('tuning', parents=[argparser], add_help=False)
    tuning_parser.set_defaults(cmd='tuning')
    tuning_parser.add_argument('--n_gcints', type=int, default=6, help='Number of GC content intervals between 0 and 1')
    # Preprocessing
    preprocessing_parser = subparsers.add_parser('preprocessing', parents=[argparser], add_help=False)
    preprocessing_parser.set_defaults(cmd='preprocessing')
    preprocessing_parser.add_argument('--db_file', type=str, help='Plasmids genes database FASTA file')    
    preprocessing_parser.add_argument('--gc_intervals', type=str, help='GC content intervals file')
    preprocessing_parser.add_argument('--out_file', type=str, help='Path to augmented dataset file')        

    return argparser.parse_args()

def main(args):    
    samples_df = _read_input(args.cmd, args.input_file)
    create_directory([args.out_dir,args.tmp_dir])
    _clean_output_files(args.cmd, samples_df, args.out_dir, args.tmp_dir)
    
    if args.cmd == 'pls_genes_db':
        create_pls_genes_db(
            args.out_dir, args.tmp_dir, samples_df
        )
    elif args.cmd == 'map_genes_to_ctgs':
        check_file(args.db_file)
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        map_pls_genes_to_contigs(
            args.out_dir, args.tmp_dir, samples_df,
            args.db_file
        )
        if args.out_file:
            _write_samples_df(samples_df, args.out_file)
            log_file(args.out_file)
    elif args.cmd == 'ground_truth':
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        create_ground_truth_files(
            args.out_dir, args.tmp_dir, samples_df
        )
        if args.out_file:
            _write_samples_df(samples_df, args.out_file)
            log_file(args.out_file)
    elif args.cmd == 'gene_density':
        create_gene_density_files(
            args.out_dir, samples_df
        )
        if args.out_file:
            _write_samples_df(samples_df, args.out_file)
            log_file(args.out_file)
    elif args.cmd == 'seeds':
        check_file(args.db_file)
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        create_seeds_parameters_file(
            args.out_dir, args.tmp_dir, samples_df,
            args.db_file
        )
    elif args.cmd == 'gc_intervals':
        create_GC_content_intervals_file(
            args.out_dir, args.tmp_dir,
            samples_df, args.n_gcints
        )
    elif args.cmd == 'gc_probabilities':
        check_file(args.gc_intervals)
        create_GC_content_probabilities_files(
            args.out_dir, args.tmp_dir,
            args.gc_intervals, samples_df
        )
        if args.out_file:
            _write_samples_df(samples_df, args.out_file)
            log_file(args.out_file)
    elif args.cmd == 'tuning':
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
            args.out_dir, args.tmp_dir,
            samples_df, args.n_gcints
        )
    elif args.cmd == 'preprocessing':
        check_file(args.db_file)
        check_file(args.gc_intervals)
        create_tmp_data_files(
            args.tmp_dir, samples_df
        )
        map_pls_genes_to_contigs(
            args.out_dir, args.tmp_dir, samples_df,
            args.db_file
        )
        create_GC_content_probabilities_files(
            args.out_dir, args.tmp_dir,
            args.gc_intervals, samples_df
        )
        create_gene_density_files(
            args.out_dir, samples_df
        )
        _write_samples_df(samples_df, args.out_file)
        log_file(args.out_file)
        
if __name__ == "__main__":
    
    args = _read_arguments()
    
    logging.basicConfig(
        filename=args.log_file,
        filemode='w',
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s'
    )
    
    main(args)
