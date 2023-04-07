import os

import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

EFAECIUM = [15,16]
ECOLI = [18,19,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,52]
KPNEUMONIAE = [62,63,64,65,66,76,85,86,87]
OTHER = [1,5,55,56,102,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,129,133]

GPLAS_IDS = EFAECIUM + ECOLI + KPNEUMONIAE
ALL_IDS = GPLAS_IDS + OTHER


thresholds = [0, 100, 1000]
eval_types = ['basepair', 'contig']

tools = {'plasbin_flow_1_1_1': 'PB_Flow', \
	 		'plasbin': 'PB', \
			'hyasp': 'HyAsP',\
			'plasmidspades': 'pSPAdes',\
			'mob_recon': 'MOB',\
			'gplas': 'gplas'}

OUT_DIR = '../eval/'

for et in eval_types:
	for th in thresholds:  

		set_1_file = 'set_1_'+ et + '_eval_gt' + str(th) + '.csv'
		set_2_file = 'set_2_'+ et + '_eval_gt' + str(th) + '.csv'

		SET_1_DICT = {}
		SET_2_DICT = {}

		sp_array_1 = []
		sp_array_2 = []

		for i in ALL_IDS:   
			if i in set(EFAECIUM):
				sp_array_1.append('E. faecium')
			elif i in set(ECOLI):
				sp_array_1.append('E. coli')
			elif i in set(KPNEUMONIAE):
				sp_array_1.append('K. pneumoniae')
			elif i == 1:
				sp_array_2.append('A. veronii')
			elif i == 5:
				sp_array_2.append('C. freundii')
			elif i == 55:
				sp_array_2.append('K. aerogenes')
			elif i == 56:
				sp_array_2.append('K. oxytoca')
			else:
				sp_array_2.append('S. enterica')

			SID = str(i)
			EVAL_DIR = '../../../evaluations/' + et + '_level/min_' + str(th) + '/sample_' + SID + '/'

			

			if i not in GPLAS_IDS:
				tools_new = list(tools.keys())[:-1]
				STATS_DICT = SET_2_DICT
			else:
				tools_new = list(tools.keys())
				STATS_DICT = SET_1_DICT

			STATS_DICT['sample_'+SID] = {}
			for tool in tools_new:
				in_file = EVAL_DIR + tool + '.txt'
				str_list = read_file(in_file)
				prec, rec, f1 = float(str_list[-3].split('\t')[1]), \
								float(str_list[-2].split('\t')[1]), \
								float(str_list[-1].split('\t')[1])
				STATS_DICT['sample_'+SID][tool + '_Prec'] = prec
				STATS_DICT['sample_'+SID][tool + '_Rec'] = rec
				STATS_DICT['sample_'+SID][tool + '_F1'] = f1
			
		SET_1_DF = pd.DataFrame.from_dict(SET_1_DICT).T
		SET_2_DF = pd.DataFrame.from_dict(SET_2_DICT).T
		SET_1_DF['Species'] = sp_array_1
		SET_2_DF['Species'] = sp_array_2
		SET_1_DF.to_csv(OUT_DIR + set_1_file)		
		SET_2_DF.to_csv(OUT_DIR + set_2_file)			











