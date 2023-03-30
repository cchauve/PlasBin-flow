#!/usr/bin/env python

# - create: create a gene database from a collection of plasmids
# - map: map a collection of genes to the contigs of an assembly graph
# - filter: remove low-quality hits from a gene-contig mapping
#
# Requirements (recursively):
#  - makeblastdb (makeblastdb / --makeblastdb; tested with BLAST+ v2.6.0)
#  - blastn (blastn / --blastn; tested with BLAST+ v2.6.0)
#  - standard UNIX tools (curl, rm)


import create_db as cd, map_genes as mg


def main():
    import argparse
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers(help = 'task to be performed', dest = 'command')
    subparsers.required = True

    # database creation
    db_parser = subparsers.add_parser('create', help = 'create gene database from collection of plasmids')
    db_parser.add_argument('genes_file', help = 'output file of genes (FASTA)')
    src = db_parser.add_mutually_exclusive_group(required = True)
    src.add_argument('-a', '--from_accession', default = cd.DEF_FROM_ACCESSION, help = 'file with one plasmid accession number per line')
    src.add_argument('-gb', '--from_genbank', default = cd.DEF_FROM_GENBANK, help = 'file with one path to plasmid GenBank file per line')
    src.add_argument('-p', '--from_plasmid_table', default = cd.DEF_FROM_PLASMID_TABLE, help = 'plasmid table downloaded from https://www.ncbi.nlm.nih.gov/genome/browse#!/plasmids/ (with all possible columns)')
    db_parser.add_argument('-k', '--keep_plasmids', default = cd.DEF_KEEP_PLASMIDS, help = 'stores the plasmid sequences (in FASTA format) if a file is specified')
    db_parser.add_argument('-d', '--dereplicate', action = 'store_true', help = 'remove duplicate genes from gene collection')
    db_parser.add_argument('-c', '--from_command_line', action = 'store_true', help = 'instead of reading the accessions (file paths), a comma-separated list of accessions (file paths) is expected after -a (-g)')
    db_parser.add_argument('-e', '--extend', action = 'store_true', help = 'add genes (and plasmids) to existing database')
    db_parser.add_argument('-r', '--released_before', default = cd.DEF_RELEASED_BEFORE, help = 'extract only genes of plasmids released before the specified data (format: YYYY-MM-DDTHH:MM:SSZ, T and Z are not replaced by digits; affects only -p)')
    db_parser.add_argument('-t', '--type', default = cd.DEF_TYPE, help = 'specify whether only RefSeq, only GenBank or both databases should be used (affects only -p)')
    db_parser.add_argument('-b', '--blacklist', default = cd.DEF_BLACKLIST, help = 'path to file containing accessions of plasmids to exclude from database (one accession per row)')
    db_parser.add_argument('-l', '--min_length', type = float, default = cd.DEF_MIN_LENGTH, help = 'minimum length of plasmids to be considered for the databases (does not apply to previously added plasmids)')
    db_parser.add_argument('-L', '--max_length', type = float, default = cd.DEF_MAX_LENGTH, help = 'maximum length of plasmids to be considered for the databases (does not apply to previously added plasmids')
    db_parser.add_argument('-m', '--min_gene_length', type = float, default = cd.DEF_MIN_GENE_LENGTH, help = 'minimum length of genes to be considered for the databases')
    db_parser.add_argument('-n', '--num_attempts', type = float, default = cd.DEF_NUM_ATTEMPTS, help = 'maximum number of attempts to properly download a GenBank file from NCBI')
    db_parser.add_argument('-v', '--verbose', action = 'store_true', help = 'print more information')

    # map genes
    map_parser = subparsers.add_parser('map', help = 'map genes to contigs')
    map_parser.add_argument('genes_file', help = 'FASTA file of plasmid genes to be mapped')
    map_parser.add_argument('mapping_file', help = 'output gene-contig mapping')
    contigs_src = map_parser.add_mutually_exclusive_group(required = True)
    contigs_src.add_argument('-f', '--from_fasta', default = mg.DEF_FROM_FASTA, help = 'file containing the contigs (in FASTA format) to which the genes should be matched')
    contigs_src.add_argument('-ag', '--from_gfa', default = mg.DEF_FROM_GFA, help = 'file containing the contigs (as part of assembly graph in GFA format) to which the genes should be matched')
    map_parser.add_argument('-c', '--clean', action = 'store_true', help = 'remove temporary files when mapping is done')
    map_parser.add_argument('-v', '--verbose', action = 'store_true', help = 'print more information')
    map_parser.add_argument('--makeblastdb', default = mg.DEF_MAKEBLASTDB_PATH, help = 'path to makeblastdb executable')
    map_parser.add_argument('--blastn', default = mg.DEF_BLASTN_PATH, help = 'path to blastn executable')

    # filter mapping
    filter_parser = subparsers.add_parser('filter', help = 'filter gene-contig mapping')
    filter_parser.add_argument('genes_file', help = 'FASTA file of plasmid genes used to create the mapping')
    filter_parser.add_argument('mapping', help = 'gene-contig mapping (BLAST output (default format 6) to be filtered')
    filter_parser.add_argument('filtered_mapping', help = 'output file containing the filtered gene-contig mapping')
    filter_parser.add_argument('-i', '--identity_threshold', type = float, default = mg.DEF_IDENTITY_THRESHOLD, help = 'minimum identity in a hit to keep it')
    filter_parser.add_argument('-l', '--length_threshold', type = float, default = mg.DEF_LENGTH_THRESHOLD, help = 'minimum fraction of query that has to be matched to keep a hit')
    filter_parser.add_argument('-f', '--find_fragmented', action = 'store_true',help = 'if set, search for fragmented hits, i.e. several short high-quality hits that together satisfy the length threshold')
    filter_parser.add_argument('-v', '--verbose', action = 'store_true', help = 'print more information')

    args = argparser.parse_args()

    if args.command == 'create':
        cd.create(args.genes_file, from_accession = args.from_accession, from_genbank = args.from_genbank,
                  from_plasmid_table = args.from_plasmid_table, keep_plasmids = args.keep_plasmids,
                  dereplicate = args.dereplicate, from_command_line = args.from_command_line, extend = args.extend,
                  released_before = args.released_before, type = args.type, blacklist = args.blacklist,
                  min_length = args.min_length, max_length = args.max_length, min_gene_length = args.min_gene_length,
                  num_attempts = args.num_attempts, verbose = args.verbose)

    if args.command == 'map':
        mg.map(args.mapping_file, args.genes_file, args.from_fasta, args.from_gfa, args.clean, args.verbose, args.makeblastdb, args.blastn)
        

    if args.command == 'filter':
        mg.filter_blast(args.genes_file, args.mapping, args.filtered_mapping, identity_threshold = args.identity_threshold,
                        length_threshold = args.length_threshold, find_fragmented = args.find_fragmented, verbose = args.verbose)


if __name__ == '__main__':
    main()