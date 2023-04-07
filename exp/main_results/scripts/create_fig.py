from sys import argv
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random

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

FIG_DIR = '../figs/'

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-sset", help="Sample set (1 or 2)")
	parser.add_argument("-et", help="Evaluation type (contig or basepair)")
	parser.add_argument("-th", help="Length threshold (0, 100 or 1000)")
	
	args = parser.parse_args()

	et = args.et
	th = args.th
	sset = args.sset

	eval_file = '../eval/set_' + str(sset) +'_' + et + '_eval_gt' + str(th) + '.csv'
	STATS = pd.read_csv(eval_file)
	STATS.rename(columns={"Unnamed: 0": "Sample"}, inplace=True)
	STATS.set_index("Sample")
	#print(STATS)

	precs, recs, f1s = [], [], []
	labels = []
	for col in STATS.columns:
		if 'Prec' in col:
			precs.append(col)
			tool = col.replace('_Prec','')
			labels.append(tools[tool])
		elif 'Rec' in col:
			recs.append(col)
		elif 'F1' in col:
			f1s.append(col)

	colors = ['red', 'orange', 'yellow', 'green', 'blue','indigo']
	tics = [1,2,3,4,5,6]
	if sset == '2':
		tics = [1,2,3,4,5]

	fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols=3, figsize=(65, 25))
    
	bplot1 = ax1.violinplot(STATS[precs], showmeans=True)
	#plot_name = 'Set '+str(sset)+' '+et+' level precision (threshold = '+str(th)+')'
	plot_name = ('Precision')
	ax1.set_title(plot_name,fontsize=30)
	for j in range(0, len(labels)):
		inds = [j + 1 + random.uniform(-0.2, 0.2) for d in STATS[precs[j]]]
		ax1.scatter(inds, STATS[precs[j]], color = colors[j], alpha = 1, s = 5, zorder = 3, label = labels[j])

	bplot2 = ax2.violinplot(STATS[recs], showmeans=True)
	#plot_name = 'Set '+str(sset)+' '+et+' level recall (threshold = '+str(th)+')'
	plot_name = ('Recall')
	ax2.set_title(plot_name,fontsize=30)
	for j in range(0, len(labels)):
		inds = [j + 1 + random.uniform(-0.2, 0.2) for d in STATS[recs[j]]]
		ax2.scatter(inds, STATS[precs[j]], color = colors[j], alpha = 1, s = 5, zorder = 3, label = labels[j])

	bplot3 = ax3.violinplot(STATS[f1s], showmeans=True)
	#plot_name = 'Set '+str(sset)+' '+et+' level F1 (threshold = '+str(th)+')'
	plot_name = ('F1')
	ax3.set_title(plot_name,fontsize=30)
	for j in range(0, len(labels)):
		inds = [j + 1 + random.uniform(-0.2, 0.2) for d in STATS[f1s[j]]]
		ax3.scatter(inds, STATS[precs[j]], color = colors[j], alpha = 1, s = 5, zorder = 3, label = labels[j])
		
	for bplot in (bplot1, bplot2, bplot3):
		for pc, color in zip(bplot['bodies'], colors):
			pc.set_facecolor(color)
			pc.set_edgecolor('black')
	
	for ax in [ax1, ax2, ax3]:
		ax.tick_params(axis='x', labelrotation = 0)
		ax.tick_params(axis='both', which='major', labelsize=20)
		ax.set_xticks(tics)
		ax.set_xticklabels(labels,fontsize=8)
	#ax3.legend(loc=7)
	
	plt.show()

	fig.savefig(FIG_DIR+'Set_'+str(sset)+'_'+et+'_level_gt'+str(th)+'.png')

	plt.close()	
