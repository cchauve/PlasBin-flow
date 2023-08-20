#!/usr/bin/env python

"""
Tuning/preprocessing for PlasBin-flow
Usage:

Creating a database of plasmid genes from the samples plasmids
python plasbin_utils.py pls_genes_db --input_file input_file --out_dir out_dir --tmp_dir tmp_dir
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  pls_fasta: gzipped FASTA plasmids file
  where in FASTA files, the header of entries are GenBank accession
- output_directory: directory where the plasmid genes database written:
  pls.genes.fasta
- tmp_directory: temporary directory, not deleted

Mapping plasmid genes to samples contigs 
python plasbin_utils.py map_genes_to_ctgs --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --db_file pls_db_file
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- output_directory: directory where the mapping files are written
  <sample>.genes_mappings.txt
- tmp_directory: temporary directory, not deleted
- pls_db_file: path to plasmid genes database file

Computing ground truth for samples
python plasbin_utils.py ground_truth --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --pid_threshold p --cov_threshold c
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  gfa: gzipped GFA file
  pls_fasta: gzipped FASTA plasmids file
- output_directory: directory where the ground truth files are written
  <sample>.ground_truth.tsv
- tmp_directory: temporary directory, not deleted
- p: percent identity threshold to define a mapping to a plasmid (default=0.95)
- c: coverage threshold to accept a blast hit  (default=0.8)

Computing GC content probabilities for samples
python plasbin_utils.py gc_probabilities --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --gc_intervals gc_intervals_file
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- output_directory: directory where the GC probabilities files are written
  <sample>.gc.tsv
- tmp_directory: temporary directory, not deleted
- gc_intervals_file: GC intervals file

Computing GC intervals from tuning samples
python plasbin_utils.py gc_intervals --input_file input_file --out_dir out_dir --tmp_dir tmp_dir
- input_file: CSV file with one line per sample and 3 required fields:
  sample: sample name
  chr_fasta: gzipped chromosome FASTA file
  pls_fasta: gzipped FASTA plasmids file
- output_directory: directory where the GC intervals files are written
  gc.txt, gc.png
- tmp_directory: temporary directory, not deleted

Computing seeds parameters
python plasbin_utils.py seeds --input_file input_file --out_dir out_dir --tmp_dir tmp_dir --db_file pls_db_file --gt_dir gt_dir
- input_file: CSV file with one line per sample and 2 required fields:
  sample: sample name
  gfa: gzipped GFA file
- output_directory: directory where the seeds parameters file is written
  seeds.txt
- tmp_directory: temporary directory, not deleted
- pls_db_file: path to plasmid genes database file
- gt_dir: directory where ground truth files <sample>.ground_truth.tsv are located

python plasbin_tuning.py tuning --input_file input_file --out_dir output_directory --tmp_dir tmp_directory 
- input_file: CSV file with one line per sample and 4 required fields:
  sample: sample name
  gfa: gzipped GFA file
  chr_fasta: gzipped chromosome FASTA file
  pls_fasta: gzipped FASTA plasmids file
  where in FASTA files, the header of entries are GenBank accession
- output_directory: directory where the tuning files are written:
  pls.genes.fasta: plasmids genes database
  gc.txt: GC content per sample
  gc.png: GC content violin plot
  seeds.txt: seed parameters
- tmp_directory: temporary directory, not deleted
"""

"""
TODO: actually generate the GC intervals file
TODO: proper logging
TODO: exceptions
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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Reading input file

def read_samples(in_csv_file):
    """
    Reads the samples information
    
    Args:
        in_csv_file (str): path to the CSV file containing tuning samples information

    Returns: 
        (DataFrame): dataframe indexed by sample id
    """
    samples_df = pd.read_csv(
        in_csv_file,
        sep = ',',
        header = 0,
        index_col = 'sample'
    )
    return samples_df

# Sample data access functions
def _get_gfa(samples_df, sample):
    return samples_df.at[sample,'gfa']
def _get_chr_fasta(samples_df, sample):
    return samples_df.at[sample,'chr_fasta']
def _get_pls_fasta(samples_df, sample):
    return samples_df.at[sample,'pls_fasta']
def _get_ground_truth(samples_df, sample):
    return samples_df.at[sample,'ground_truth']
def _set_ground_truth(samples_df, sample, gt_file):
    samples_df.at[sample,'ground_truth'] = gt_file

# Auxiliary functions

def gunzip_fasta(in_fasta_file, out_fasta_file):
    records = []
    with gzip.open(in_fasta_file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
    with open(out_fasta_file, 'w') as out_file:
        SeqIO.write(records, out_file, 'fasta')
        
def gfa2fasta(in_gfa_file, out_fasta_file):
    """
    Convert a GFA file into a FASTA file

    Args:
       in_gfa_file (str): path to input gzipped GFA file
       out_fasta_file (str): path to output FASTA file
    """
    gfa_ctg_seqs = {}
    with open(in_gfa_file, 'r') as in_file, open(out_fasta_file, 'w') as out_file:
        for line in in_file.readlines():
            line_split = line.strip().split('\t')
            if line[0] == 'S':
                ctg_id,ctg_data = line_split[1],line_split[2:]
                ctg_seq = ctg_data[0]
                out_file.write(f'>{ctg_id}\n{ctg_seq}\n')	   

def _clean_files(files2clean):
    for in_file in files2clean:
        if os.path.isfile(in_file):
            os.remove(in_file)
        
def _create_directory(in_dir_list):
    for in_dir in in_dir_list:
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

def _run_cmd(cmd, verbose):
    process = subprocess.run(cmd, capture_output=True, text=True)
    if verbose:
        print(f'#STDOUT:{process.stdout}#')
        print(f'#STDERR:{process.stderr}#')

def _read_seq_len(in_fasta_file):
    seq_len = {}
    with open(in_fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seq_len[record.id] = len(record.seq)
    return seq_len

# Temporary file names

TMP_PLS_GB_FILE='pls.genbank.txt'
def tmp_pls_gb_file(tmp_dir):
    return os.path.join(tmp_dir, TMP_PLS_GB_FILE)
def tmp_fastas_file(tmp_dir, file_type):
    return {
        'chr': os.path.join(tmp_dir, 'chr.fasta.txt'),
        'pls': os.path.join(tmp_dir, 'pls.fasta.txt')
    }[file_type]
def tmp_gfa_file(tmp_dir, sample):
    return os.path.join(tmp_dir, f'{sample}.gfa')
def tmp_gfa_fasta_file(tmp_dir, sample):
    return os.path.join(tmp_dir, f'{sample}.gfa.fasta')
def tmp_pls_fasta_file(tmp_dir, sample):
    return os.path.join(tmp_dir, f'{sample}.pls.fasta')
def tmp_pls_blastdb_file(tmp_dir, sample):
    return os.path.join(tmp_dir, f'{sample}.pls.fasta.db')
def tmp_pls_mappings_file(tmp_dir, sample):
    return os.path.join(tmp_dir, f'{sample}.pls_mappings.txt')
def tmp_seeds_input_file(tmp_dir):
    return os.path.join(tmp_dir, 'seeds_input.csv')

# Output file names

PLS_GENES_DB_FILE='pls.genes.fasta'
def pls_genes_db_file(out_dir, out_file=PLS_GENES_DB_FILE):
    return os.path.join(out_dir, out_file)
def ground_truth_file(out_dir, sample):
    return os.path.join(out_dir, f'{sample}.ground_truth.tsv')
SEEDS_PARAMETERS_FILE='seeds.txt'
def seeds_parameters_file(out_dir, out_file=SEEDS_PARAMETERS_FILE):
    return os.path.join(out_dir, out_file)

GC_FILE_PREFIX='gc'
def gc_txt_file(out_dir, out_file=f'{GC_FILE_PREFIX}.txt'):
    return os.path.join(out_dir, out_file)
def gc_png_file(out_dir, out_file=f'{GC_FILE_PREFIX}.png'):
    return os.path.join(out_dir, out_file)
def gc_proba_file(out_dir, sample):
    return os.path.join(out_dir, f'{sample}.gc.tsv')

def genes_mappings_file(tmp_dir, sample):
    return os.path.join(tmp_dir, f'{sample}.genes_mappings.txt')

def create_tmp_unzipped_gfa_fasta_files(tmp_dir, samples_df, verbose=True):
    """
    Creates a ungzipped GFA and FASTA file for each sample

    Args:
       samples_df (DataFrame): samples dataframe
       tmp_dir (str): path to temporary directory

    Returns:
       None, creates files tmp_gfa_file(tmp_dir, sample) and tmp_gfa_fasta_file(tmp_dir, sample) for each sample
    """
    for sample in samples_df.index:
        if verbose: print(f'ACTION\tcopy assembly files for {sample}')
        gfagz_file = _get_gfa(samples_df,sample)
        gfa_file = tmp_gfa_file(tmp_dir, sample)
        if not os.path.isfile(gfa_file):
            with gzip.open(gfagz_file) as in_file:
                with open(gfa_file, 'wb') as out_file:
                    shutil.copyfileobj(in_file, out_file)
        if verbose: print(f'FILE\t{tmp_gfa_file(tmp_dir, sample)}')
        gfa_fasta_file = tmp_gfa_fasta_file(tmp_dir, sample)
        if not os.path.isfile(gfa_fasta_file):
            gfa2fasta(
                tmp_gfa_file(tmp_dir, sample),
                gfa_fasta_file
            )
        if verbose: print(f'FILE\t{tmp_gfa_fasta_file(tmp_dir, sample)}')


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
      None, creates the file ground_truth_file(out_dir, sample)
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
        tmp_pls_mappings_file(tmp_dir, sample),
        sep = '\t', names = col_names, dtype = str
    )
        
    ctg_len = _read_seq_len(tmp_gfa_fasta_file(tmp_dir, sample))        
    pls_len = _read_seq_len(tmp_pls_fasta_file(tmp_dir, sample))
    
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
    with open(ground_truth_file(out_dir, sample), "w") as out_file:
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
        out_dir, tmp_dir, samples_df, pid_threshold=0.95, cov_threshold=0.8,
        verbose=True
):
    """
    Creates a TSV ground truth file per sample

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe

    Returns:
       None, creates files ground_truth_file(out_dir, sample) for each sample  
    """
    for sample in samples_df.index:
        if verbose: print(f'ACTION\tground truth for {sample}')
        pls_fastagz_file = _get_pls_fasta(samples_df, sample)
        pls_fasta_file = tmp_pls_fasta_file(tmp_dir, sample)
        pls_blastdb = tmp_pls_blastdb_file(tmp_dir, sample)
        gfa_fasta_file = tmp_gfa_fasta_file(tmp_dir, sample)
        mappings_file = tmp_pls_mappings_file(tmp_dir, sample)
        gunzip_fasta(pls_fastagz_file, pls_fasta_file)
        if verbose: print(f'ACTION\tcompute blast database for {pls_fasta_file}')
        cmd1 = [
            'makeblastdb',
            '-in', pls_fasta_file,
            '-dbtype', 'nucl',
            '-out', pls_blastdb
        ]
        _run_cmd(cmd1, verbose)
        if verbose: print(f'FILE\t{pls_blastdb}')        
        if verbose: print(f'ACTION\tmap {gfa_fasta_file} to {pls_blastdb}')        
        cmd2 = [
            'blastn', '-task', 'megablast',
            '-query', gfa_fasta_file,
            '-db', pls_blastdb,
            '-out', mappings_file,
            '-outfmt', '6'
        ]
        _run_cmd(cmd2, verbose)
        if verbose: print(f'FILE\t{mappings_file}')
        if verbose: print(f'ACTION\tcompute ground truth file')                
        _compute_ground_truth(
            out_dir, tmp_dir, sample, pid_threshold, cov_threshold
        )
        _set_ground_truth(
            samples_df, sample, ground_truth_file(out_dir, sample)
        )
        if verbose: print(f'FILE\t{ground_truth_file(out_dir, sample)}')
    
def create_pls_genes_db(
        out_dir, tmp_dir, samples_df, out_file_name=PLS_GENES_DB_FILE,
        verbose=True
):
    """
    Creates a plasmid genes database

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       out_file_name (str): name of output file

    Returns:
       None, creates file pls_genes_db_file(out_dir, out_file_name)
    """
    if verbose: print(f'ACTION\tcompute plasmid genes database')
    with open(tmp_pls_gb_file(tmp_dir), 'w') as out_file:
        for sample in samples_df.index:
            if verbose: print(f'ACTION\trecord plasmid GenBank accession for {sample}')
            with gzip.open(_get_pls_fasta(samples_df, sample), 'rt') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    out_file.write(f'{record.id}\n')
    if verbose: print(f'FILE\t{tmp_pls_gb_file(tmp_dir)}')
    if verbose: print(f'ACTION\tprocess {tmp_pls_gb_file(tmp_dir)}')
    cmd = [
        'python', 'get_gd.py', 'create',
        pls_genes_db_file(out_dir, out_file_name),
        '-a', tmp_pls_gb_file(tmp_dir)
    ]
    _run_cmd(cmd, verbose)
    if verbose: print(f'FILE\t{pls_genes_db_file(out_dir, out_file_name)}')

def map_pls_genes_to_contigs(out_dir, tmp_dir, samples_df, db_file, verbose=True):
    """
    Map plasmid genes to samples contigs

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       db_file (str): path to plasmid genes database file

    Returns:
       None, creates files genes_mappings_file(out_dir, sample)nfor each sample
    """
    for sample in samples_df.index:
        mappings_file = genes_mappings_file(out_dir, sample)
        if verbose: print(f'ACTION\tmapping {sample} to {db_file}')
        cmd = [
            'python', 'get_gd.py', 'map',
            db_file,
            mappings_file,
            '-a', tmp_gfa_file(tmp_dir, sample)
        ]
        _run_cmd(cmd, verbose)
        if verbose: print(f'FILE\t{mappings_file}')
        
def create_GC_content_intervals_file(
        out_dir, tmp_dir, samples_df, out_file_prefix=GC_FILE_PREFIX,
        verbose=True
):
    """
    Creates GC content intervals files

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       out_file_prefix (str): prefix of output file names

    Returns:
       None, creates files gc_txt_file(out_dir) and gc_png_file(out_dir)
    """
    if verbose: print(f'ACTION\tcompute GC content intervals files')
    file_access_fun = {
        'chr': _get_chr_fasta,
        'pls': _get_pls_fasta
    }
    for file_type in ['chr','pls']:
        with open(tmp_fastas_file(tmp_dir, file_type), 'w') as out_file:
            if verbose: print(f'ACTION\trecord {file_type} files paths')
            for sample in samples_df.index:
                fasta_file = file_access_fun[file_type](samples_df, sample)
                out_file.write(f'{fasta_file}\n')
            if verbose: print(f'FILE\t{tmp_fastas_file(tmp_dir, file_type)}')
    out_txt_file = gc_txt_file(out_dir, out_file=f'{out_file_prefix}.txt')
    out_png_file = gc_png_file(out_dir, out_file=f'{out_file_prefix}.png')
    if verbose: print(f'ACTION\tcompute output files')
    cmd = [
        'python', 'analyse_GC_content.py',
        '--chr', tmp_fastas_file(tmp_dir, 'chr'),
        '--pls', tmp_fastas_file(tmp_dir, 'pls'),        
        '--out', out_txt_file,
        '--vplot', out_png_file
    ]
    _run_cmd(cmd, verbose)
    if verbose: print(f'FILE\t{out_txt_file}\nFILE\t{out_png_file}')

def create_GC_content_probabilities_file(
    out_dir, tmp_dir, gc_intervals_file, samples_df,verbose=True
): 
    """
    Creates GC content probabilities files

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe

    Returns:
       None, creates files gc_proba_file(out_dir, sample) for each sample
    """
    if verbose: print(f'ACTION\tcompute GC content probabilities files')
    for sample in samples_df.index:
        gfa_file = tmp_gfa_file(tmp_dir, sample)
        gcp_file = gc_proba_file(out_dir, sample)
        cmd = [
            'python', 'get_gc_probs.py',
            '-ag', gfa_file,
            '-outfile', gcp_file,
            '-gcint', gc_intervals_file
        ]
        _run_cmd(cmd, verbose)
        if verbose: print(f'FILE\t{gcp_file}')

def create_seeds_parameters_file(
        out_dir, tmp_dir, samples_df,
        db_file, gt_dir, out_file_name=SEEDS_PARAMETERS_FILE,
        verbose=True
):
    """
    Creates a file containing the parameters defining seeds

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       samples_df (DataFrame): samples dataframe
       db_file (str): path to plasmid genes database file
       gt_dir (str): directory containing groud truth files
       out_file_name (str): name of output file

    Returns:
       None, creates file seeds_parameters_file(out_dir, out_file_name)
    """
    if verbose: print(f'ACTION\tcompute seeds parameters')    
    map_pls_genes_to_contigs(
        tmp_dir, tmp_dir, samples_df, db_file, verbose=verbose
    )
    with open(tmp_seeds_input_file(tmp_dir), 'w') as out_file:
        for sample in samples_df.index:
            mappings_file = genes_mappings_file(tmp_dir, sample)
            gt_file = ground_truth_file(gt_dir, sample)
            gfa_fasta_file = tmp_gfa_fasta_file(tmp_dir, sample)
            out_file.write(
                f'{sample},{gfa_fasta_file},{mappings_file},{gt_file}\n'
            )
    if verbose: print(f'FILE\t{tmp_seeds_input_file(tmp_dir)}')
    cmd = [
        'python', 'analyse_seed_eligibility.py',
        '--paths', tmp_seeds_input_file(tmp_dir),
        '--out', seeds_parameters_file(out_dir, out_file_name)
    ]
    _run_cmd(cmd, verbose)
    if verbose: print(f'FILE\t{seeds_parameters_file(out_dir, out_file_name)}')


def main():
    argparser = argparse.ArgumentParser(description='PlasBin-flow utils')
    argparser.add_argument('--verbose', type=bool, default=True, help='Verbose mode')
    argparser.add_argument('--input_file', type=str, help='Samples CSV file')
    argparser.add_argument('--out_dir', type=str, help='Output directory')    
    argparser.add_argument('--tmp_dir', type=str, help='Temporary directory')
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
    gt_parser.add_argument('--pid_threshold', type=float, default=0.95, help='Percent identity threshold in [0,1]')
    gt_parser.add_argument('--cov_threshold', type=float, default=0.8, help='Percent coverage threshold in [0,1]')    
    # Computing seeds parameters
    seeds_parser = subparsers.add_parser('seeds', parents=[argparser], add_help=False)
    seeds_parser.set_defaults(cmd='seeds')
    seeds_parser.add_argument('--db_file', type=str, help='Plasmids genes database FASTA file')
    seeds_parser.add_argument('--gt_dir', type=str, help='Directory containing ground truth files')    
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

    if args.cmd == 'pls_genes_db':
        samples_df = read_samples(args.input_file)
        files2clean = [
            pls_genes_db_file(args.out_dir)
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_pls_genes_db(
            args.out_dir, args.tmp_dir, samples_df,
            verbose=args.verbose
        )

    elif args.cmd == 'map_genes_to_ctgs':
        samples_df = read_samples(args.input_file)
        files2clean = [
            genes_mappings_file(args.out_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_unzipped_gfa_fasta_files(
            args.tmp_dir, samples_df, verbose=args.verbose
        )
        map_pls_genes_to_contigs(
            args.out_dir, args.tmp_dir, samples_df, args.db_file,
            verbose=args.verbose
        )

    elif args.cmd == 'ground_truth':
        samples_df = read_samples(args.input_file)
        files2clean = [
            ground_truth_file(args.out_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_unzipped_gfa_fasta_files(
            args.tmp_dir, samples_df, verbose=args.verbose
        )
        create_ground_truth_files(
            args.out_dir, args.tmp_dir, samples_df,
            verbose=args.verbose
        )

    elif args.cmd == 'seeds':
        samples_df = read_samples(args.input_file)
        files2clean = [seeds_parameters_file(args.out_dir)]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_unzipped_gfa_fasta_files(
            args.tmp_dir, samples_df, verbose=args.verbose
        )
        create_seeds_parameters_file(
            args.out_dir, args.tmp_dir, samples_df,
            args.db_file, args.gt_dir, verbose=args.verbose
        )

    elif args.cmd == 'gc_intervals':
        samples_df = read_samples(args.input_file)
        files2clean = [
            gc_txt_file(args.out_dir),
            gc_png_file(args.out_dir)
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_GC_content_intervals_file(
            args.out_dir, args.tmp_dir, samples_df,
            verbose=args.verbose
        )

    elif args.cmd == 'gc_probabilities':
        samples_df = read_samples(args.input_file)
        files2clean = [
            gc_proba_file(args.out_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_GC_content_probabilities_file(
            args.out_dir, args.tmp_dir,
            args.gc_intervals, samples_df,
            verbose=args.verbose
        )
        
    elif args.cmd == 'tuning':
        samples_df = read_samples(args.input_file)
        files2clean = [
            pls_genes_db_file(args.out_dir),
            seeds_parameters_file(args.out_dir),
            gc_txt_file(args.out_dir),
            gc_png_file(args.out_dir)
        ] + [
            ground_truth_file(args.tmp_dir, sample)
            for sample in samples_df.index
        ]
        _clean_files(files2clean)
        _create_directory([args.out_dir,args.tmp_dir])
        create_tmp_unzipped_gfa_fasta_files(
            args.tmp_dir, samples_df, verbose=args.verbose
        )
        create_pls_genes_db(
            args.out_dir, args.tmp_dir, samples_df,
            verbose=args.verbose
        )
        create_ground_truth_files(
            args.out_dir, args.tmp_dir, samples_df,
            verbose=args.verbose
        )
        create_seeds_parameters_file(
            args.out_dir, args.tmp_dir, samples_df,
            pls_genes_db_file(args.out_dir), args.out_dir,
            verbose=args.verbose
        )
        create_GC_content_intervals_file(
            args.out_dir, args.tmp_dir, samples_df,
            verbose=args.verbose
        )
        
if __name__ == "__main__":
    main()
