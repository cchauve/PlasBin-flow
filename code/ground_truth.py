import pandas as pd

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

def compute_ground_truth_file(
        out_dir, tmp_dir, sample, pls_mappings_file,
        ctg_len, pls_len, pid_threshold, cov_threshold,
        ground_truth_file):
    """
    Computes the ground truth file for a sample

    Args:
       out_dir (str): path to output directory
       tmp_dir (str): path to temporary directory
       sample (str): sample -> sample id
       pls_mappings_file (str): path to file of mappings of contigs to true plasmids
       ctg_len (Dictionary): contig id -> contig length
       pls_len (Dictionary): plasmid id -> plasmid length
       pid_threshold (float): percent identity threshold
       cov_threshold (float): coverage threshold
       ground_truth_file (str): path to the ground truth file to write

    Returns:
      None, creates the file _ground_truth_file(out_dir, sample)
    """
    col_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore"
    ]  # outfmt 6
    hits = pd.read_csv(
        pls_mappings_file,
        sep = '\t', names = col_names, dtype = str
    )
        
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
    with open(ground_truth_file, "w") as out_file:
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
