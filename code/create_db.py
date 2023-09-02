#!/usr/bin/env python

# Modified from the code for the command 'create' of HyAsP.
# https://github.com/cchauve/HyAsP
#
# Takes a list of plasmids with annotated genes and extracts all genes into a FASTA file.
# The plasmids can be provided as GenBank accession numbers (in this case,
# the GenBank files are temporarily downloaded).
#
# Additional filtering options are available for this source:
#  - blacklisting of plasmids not to be included in the database
#
# There are further options allowing to dereplicate the genes (based on their sequence),
# to consider only plasmids in a certain length range
#
# Requirements:
# standard UNIX tools (curl, rm)


import math
import os
import pandas as pd
import logging

from Bio import Seq, SeqIO

from log_errors_utils import (
    process_exception,
    run_cmd_redirect,
    run_cmd
)

# default values / constants
DEF_DEREPLICATE = True
DEF_BLACKLIST = ''
DEF_MIN_LENGTH = 0
DEF_MAX_LENGTH = math.inf
DEF_MIN_GENE_LENGTH = 0
DEF_NUM_ATTEMPTS = 25

# reads GenBank file and extracts all genes from it
def extract_all_genes(gb_file):
    gene_collection = []
    for gb_record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):
        # find all genes
        genes = []
        for feature in gb_record.features:
            if feature.type == 'gene':
                genes.append(feature)
        # create FASTA entry per gene
        i = 0
        for g in genes:
            if 'locus_tag' in g.qualifiers:
                gene_collection.append((g.qualifiers['locus_tag'][0], str(g.extract(gb_record.seq))))
            elif 'gene' in g.qualifiers:
                gene_qualifier = g.qualifiers['gene'][0]
                gene_collection.append((f'gene{i}_{gene_qualifier}', str(g.extract(gb_record.seq))))
            else:
                logging.warning(f'No identifier for  gene {g} in {gb_file}')
            i += 1
    return gene_collection

# reads GenBank file and extracts the sequence of the plasmid
def extract_seq(gb_file):
    seqs = []
    for gb_record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):
        try:
            seqs.append((gb_record.name, str(gb_record.seq)))
        except Seq.UndefinedSequenceError as e:
            process_exception(f'Undefined sequence {gb_record.name} in {gb_file}: {e}')
    return seqs

# creates gene database (and plasmids database) from given list of plasmids
def create_db(
        genes_file, sources, blacklist=DEF_BLACKLIST, dereplicate=DEF_DEREPLICATE,
        min_length = DEF_MIN_LENGTH, max_length = DEF_MAX_LENGTH,
        min_gene_length = DEF_MIN_GENE_LENGTH,
        num_attempts = DEF_NUM_ATTEMPTS
):
    temp = f'{genes_file}_tmp_gb'

    # remove blacklisted accessions
    sources = [src for src in sources if src not in blacklist]
    logging.info(f'{len(sources)} plasmid accessions after removing blacklisted plasmids.')

    # write old (dereplicated) contents and add genes from new plasmids
    num_genes = 0
    num_length_discarded = 0
    num_attempts_discarded = 0
    gene_seqs = set()
    with open(genes_file, 'w') as out_genes:
        for i, src in enumerate(sources, start = 1):
            plasmid_seqs = []
            res = -1
            cnt_failures = 0
            logging.info(f'ACTION\tDownloading plasmid with genBank accession {src}.')
            while res != 0 and cnt_failures < num_attempts:
                cmd = [
                    'curl',
                    f'https://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={src}&rettype=gbwithparts&retmode=text'
                ]
                res = run_cmd_redirect(cmd, temp, num_attempts=1, exit_on_error=False)
                if res != 0:
                    logging.warning(
                        f'Download of plasmid {src} failed. Trying again.'
                    )
                    cnt_failures += 1
                else:
                    plasmid_seqs = extract_seq(temp)
                    if len(plasmid_seqs) == 0:
                        logging.warning(
                            f'No sequence for plasmid was found in {src}. Download seems to be faulty. Repeating download.'
                        )
                        res = -1
                        cnt_failures += 1

            if cnt_failures == num_attempts:
                logging.warning(
                    f'Plasmid {src} could not be downloaded properly (within {num_attempts} attempts). Continuing with next plasmid.'
                )
                num_attempts_discarded += 1
                plasmid_seqs = None
                genes = []
            else:
                genes = extract_all_genes(temp)

            if plasmid_seqs is not None and len(plasmid_seqs) > 0:
                if len(plasmid_seqs) > 1:
                    logging.warning(
                        f'More than one plasmid sequence in {src}. Using only the first one.'
                    )
                plasmid_name, plasmid_seq = plasmid_seqs[0]

                if min_length <= len(plasmid_seq) <= max_length:
                    for id, seq in genes:
                        if ((not dereplicate) or (seq not in gene_seqs)) and len(seq) >= min_gene_length:
                            gene_seqs.add(seq)
                            out_genes.write('>%s\n%s\n' % (id, seq))
                            num_genes += 1
                else:
                    num_length_discarded += 1

    if num_length_discarded > 0:
        logging.warning(f'{num_length_discarded} plasmids were discarded for their length.')
    if num_attempts_discarded > 0:
        logging.warning(f'{num_attempts_discarded} plasmids were discarded for reaching the download-attempt limit.')

    logging.info(f'Database comprises {num_genes} genes.')
    _ = run_cmd(['rm', '-f', temp])

# print configuration of database generation
def log_config(
        from_accession, genes_file, dereplicate, blacklist,
        min_length , max_length, min_gene_length,
        num_attempts
):
    blacklist_str = ','.join(blacklist) if len(blacklist) > 0 else '(empty)'
    logging.info(
        f'''############################################
### Configuration of database generation ###
>>> Input / output
Plasmid source: Accession
Genes file: {genes_file}
>>> Filtering
Minimum plasmid length: {min_length}
Maximum plasmid length: {max_length}
Minimum gene length: {min_gene_length}
Blacklist: {blacklist_str}
>>> Other options
Dereplicate: {dereplicate}
Maximum number of download attempts: {num_attempts}
############################################''')


# choose correct method depending on inputs and make sure that a valid selections of options is provided
def create(
        genes_file, from_accession,
        blacklist = DEF_BLACKLIST, dereplicate=DEF_DEREPLICATE,
        min_length = DEF_MIN_LENGTH, max_length = DEF_MAX_LENGTH,
        min_gene_length = DEF_MIN_GENE_LENGTH,
        num_attempts = DEF_NUM_ATTEMPTS
):

    blacklisted = []
    if blacklist != '':
        with open(blacklist, 'r') as in_file:
            for line in in_file:
                blacklisted.append(line.strip())

    log_config(from_accession, genes_file, dereplicate, blacklisted, min_length, max_length, min_gene_length, num_attempts)

    with open(from_accession, 'r') as infile:
        sources = [
            line.strip()
            for line in infile
            if line.strip()
        ]
    create_db(genes_file, sources, blacklisted, dereplicate, min_length, max_length, min_gene_length, num_attempts)
