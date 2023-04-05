import os

weights = [[1,1,2],[1,2,1],[1,2,2],[2,1,1],[2,1,2],[2,2,1]]

PBF_DIR = '../../code'
EXAMPLE_DIR = '../../example'

############################################
for wts in weights:
	#print(BACTERIUM, SAMPLE, wts)
	w1, w2, w3 = wts[0], wts[1], wts[2]
	OUT_DIR = EXAMPLE_DIR + '/' + str(w1) + '_' + str(w2) + '_' + str(w3)

	os.system('python ' + PBF_DIR + '/plasmids_flow.py -ag '+ EXAMPLE_DIR +'/input/assembly.gfa \
					-gc '+ EXAMPLE_DIR +'/input/gc_probs.tsv \
					-map '+ EXAMPLE_DIR +'/input/gene_contig_mapping.tsv \
					-alpha1 '+ str(w1) + ' -alpha2 '+ str(w2) + ' -alpha3 '+ str(w3) +\
					' -outdir '+ OUT_DIR + ' -outfile pbf_bins.out')









