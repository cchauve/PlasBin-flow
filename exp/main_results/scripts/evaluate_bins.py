#!/usr/bin/python

#USAGE: python evaluate_bins.py -gt GROUND_TRUTH_CONTIG_MAPPING -out OUT_DIR -eval_dir EVAL_DIR -tool plasbin-flow -th LENGTH_THRESHOLD

from __future__ import division
import argparse
import pandas as pd

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

def eval_results(pred_dict, pls_dict, len_dict, th_len, ue, we):
	ue.write("##### predictions #####\n")
	ue.write(">number of predicted plasmids: %i\n" % len(pred_dict))
	we.write("##### predictions #####\n")
	we.write(">number of predicted plasmids: %i\n" % len(pred_dict))

	for pred in pred_dict:
		ue.write(pred+"\n")
		we.write(pred+"\n")
		
	unwtd_recall = {}
	unwtd_precision = {}
	wtd_recall = {}
	wtd_precision = {}


	for ref_pls in pls_dict:
		unwtd_recall[ref_pls] = {}
		unwtd_recall[ref_pls]['Val'] = 0
		unwtd_recall[ref_pls]['Bin'] = ''
		unwtd_recall[ref_pls]['Common'] = 0
		unwtd_recall[ref_pls]['Total'] = 0

		wtd_recall[ref_pls] = {}
		wtd_recall[ref_pls]['Val'] = 0
		wtd_recall[ref_pls]['Bin'] = ''
		wtd_recall[ref_pls]['Common'] = 0
		wtd_recall[ref_pls]['Total'] = 0

		nref_ctgs = 0
		lenref_ctgs = 0

		for ctg in pls_dict[ref_pls]:
			if len_dict[ctg] >= th_len:
				nref_ctgs += 1
				lenref_ctgs += len_dict[ctg]
		
		for pred_pls in pred_dict:
			common_ctgs = set(pred_dict[pred_pls]).intersection(set(pls_dict[ref_pls]))
			ncommon_ctgs = 0
			lencommon_ctgs = 0
			for ctg in common_ctgs:
				if len_dict[ctg] >= th_len:
					ncommon_ctgs += 1
					lencommon_ctgs += len_dict[ctg]
			
			unwtd_rec = 0
			wtd_rec = 0
			if nref_ctgs >= 1:
				unwtd_rec = ncommon_ctgs / nref_ctgs
				wtd_rec = lencommon_ctgs / lenref_ctgs

			if unwtd_rec >= unwtd_recall[ref_pls]['Val']:
				unwtd_recall[ref_pls]['Val'] = unwtd_rec
				unwtd_recall[ref_pls]['Bin'] = pred_pls
				unwtd_recall[ref_pls]['Common'] = ncommon_ctgs
				unwtd_recall[ref_pls]['Total'] = nref_ctgs

			if wtd_rec >= wtd_recall[ref_pls]['Val']:
				wtd_recall[ref_pls]['Val'] = wtd_rec
				wtd_recall[ref_pls]['Bin'] = pred_pls
				wtd_recall[ref_pls]['Common'] = lencommon_ctgs
				wtd_recall[ref_pls]['Total'] = lenref_ctgs


	for pred_pls in pred_dict:
		unwtd_precision[pred_pls] = {}
		unwtd_precision[pred_pls]['Val'] = 0
		unwtd_precision[pred_pls]['Ref'] = ''
		unwtd_precision[pred_pls]['Common'] = 0
		unwtd_precision[pred_pls]['Total'] = 0

		wtd_precision[pred_pls] = {}
		wtd_precision[pred_pls]['Val'] = 0
		wtd_precision[pred_pls]['Ref'] = ''
		wtd_precision[pred_pls]['Common'] = 0
		wtd_precision[pred_pls]['Total'] = 0

		npred_ctgs = 0
		lenpred_ctgs = 0

		for ctg in pred_dict[pred_pls]:
			if len_dict[ctg] >= th_len:
				npred_ctgs += 1
				lenpred_ctgs += len_dict[ctg]

		for ref_pls in pls_dict:
			common_ctgs = set(pls_dict[ref_pls]).intersection(set(pred_dict[pred_pls]))
			ncommon_ctgs = 0
			lencommon_ctgs = 0
			for ctg in common_ctgs:
				if len_dict[ctg] >= th_len:
					ncommon_ctgs += 1
					lencommon_ctgs += len_dict[ctg]
			unwtd_prec = 0
			wtd_prec = 0
			if npred_ctgs >= 1:
				unwtd_prec = ncommon_ctgs / npred_ctgs
				wtd_prec = lencommon_ctgs / lenpred_ctgs

			if unwtd_prec >= unwtd_precision[pred_pls]['Val']:
				unwtd_precision[pred_pls]['Val'] = unwtd_prec
				unwtd_precision[pred_pls]['Ref'] = ref_pls
				unwtd_precision[pred_pls]['Common'] = ncommon_ctgs
				unwtd_precision[pred_pls]['Total'] = npred_ctgs

			if wtd_prec >= wtd_precision[pred_pls]['Val']:
				wtd_precision[pred_pls]['Val'] = wtd_prec
				wtd_precision[pred_pls]['Ref'] = ref_pls
				wtd_precision[pred_pls]['Common'] = lencommon_ctgs
				wtd_precision[pred_pls]['Total'] = lenpred_ctgs

	### UNWEIGHTED STATS ###

	ue.write("###Precision details###\n")
	ue.write('>Precision: Proportion of correctedly identified contigs for each prediction\n')
	ue.write("#Predicted_bin\tPrecision\tReference_plasmid\n")

	unwtd_overall_common = 0
	unwtd_overall_total = 0
	unwtd_overall_prec = 0
	for x in unwtd_precision:
		val = float("{:.4f}".format(unwtd_precision[x]['Val']))
		ref = unwtd_precision[x]['Ref']
		print(x, val, ref)
		unwtd_overall_common += unwtd_precision[x]['Common']
		unwtd_overall_total += unwtd_precision[x]['Total']
		ue.write(x + '\t' + str(val) + '\t' + ref + "\n")
	if unwtd_overall_total != 0:
		unwtd_overall_prec = unwtd_overall_common/unwtd_overall_total
	ue.write("\n")

	ue.write("###Recall details###\n")
	ue.write('>Recall: Proportion of correctedly identified contigs for each reference\n')
	ue.write("#Reference_plasmid\tRecall\tPredicted_bin\n")

	unwtd_overall_common = 0
	unwtd_overall_total = 0
	unwtd_overall_rec = 0
	for x in unwtd_recall:
		val = float("{:.4f}".format(unwtd_recall[x]['Val']))
		bin = unwtd_recall[x]['Bin']
		print(x, val, bin)
		unwtd_overall_common += unwtd_recall[x]['Common']
		unwtd_overall_total += unwtd_recall[x]['Total']
		ue.write(x + '\t' + str(val) + '\t' + bin + "\n")
	if unwtd_overall_total != 0:
		unwtd_overall_rec = unwtd_overall_common/unwtd_overall_total
	print("\n")
	ue.write("\n")

	unwtd_f1 = 0
	if (unwtd_overall_prec + unwtd_overall_rec) != 0:
		unwtd_f1 = 2*unwtd_overall_prec*unwtd_overall_rec / (unwtd_overall_prec + unwtd_overall_rec)
	
	ue.write("###Overall details###\n")
	ue.write('>Final statistics:\n')
	ue.write('Precision\t' + str(unwtd_overall_prec) + '\n')
	ue.write('Recall\t' + str(unwtd_overall_rec) + '\n')
	ue.write('F1\t' + str(unwtd_f1)  + '\n')


	### WEIGHTED STATS ###
	we.write("###Precision details###\n")
	we.write('>Precision: Proportion of correctedly identified contigs for each prediction\n')
	we.write("#Predicted_bin\tPrecision\tReference_plasmid\n")

	wtd_overall_common = 0
	wtd_overall_total = 0
	wtd_overall_prec = 0
	for x in wtd_precision:
		val = float("{:.4f}".format(wtd_precision[x]['Val']))
		ref = wtd_precision[x]['Ref']
		print(x, val, ref)
		wtd_overall_common += wtd_precision[x]['Common']
		wtd_overall_total += wtd_precision[x]['Total']
		we.write(x + '\t' + str(val) + '\t' + ref + "\n")
	if wtd_overall_total != 0:
		wtd_overall_prec = wtd_overall_common/wtd_overall_total
	we.write("\n")

	we.write("###Recall details###\n")
	we.write('>Recall: Proportion of correctedly identified contigs for each reference\n')
	we.write("#Reference_plasmid\tRecall\tPredicted_bin\n")

	wtd_overall_common = 0
	wtd_overall_total = 0
	wtd_overall_rec = 0
	for x in wtd_recall:
		val = float("{:.4f}".format(wtd_recall[x]['Val']))
		bin = wtd_recall[x]['Bin']
		print(x, val, bin)
		wtd_overall_common += wtd_recall[x]['Common']
		wtd_overall_total += wtd_recall[x]['Total']
		we.write(x + '\t' + str(val) + '\t' + bin + "\n")
	if wtd_overall_total != 0:
		wtd_overall_rec = wtd_overall_common/wtd_overall_total
	print("\n")
	we.write("\n")

	wtd_f1 = 0
	if (wtd_overall_prec + wtd_overall_rec) != 0:
		wtd_f1 = 2*wtd_overall_prec*wtd_overall_rec / (wtd_overall_prec + wtd_overall_rec)
	
	we.write("###Overall details###\n")
	we.write('>Final statistics:\n')
	we.write('Precision\t' + str(wtd_overall_prec) + '\n')
	we.write('Recall\t' + str(wtd_overall_rec) + '\n')
	we.write('F1\t' + str(wtd_f1)  + '\n')


def main():
	argparser = argparse.ArgumentParser()
	argparser.add_argument("-gt", "-ground_truth", help = "Path to ground truth file")
	argparser.add_argument("-out", "-out_dir", help = "Path to method output directory")
	argparser.add_argument("-eval_dir", help = "Path to evaluations directory")
	argparser.add_argument("-tool", help = "Method to be evaluated")
	argparser.add_argument("-th", "-len_th", default = 0, help = "Length threshold")

	args = argparser.parse_args()
	len_th = int(args.len_th)
	out = args.out_dir

	#Read ground truth file
	string_list = read_file(args.ground_truth)
	#print(len(string_list))

	pls_ids = []

	#chr_dict = {}
	pls_dict = {}

	len_dict = {}
	for line in string_list[1:]:
		line = line.split('\t')
		#print(line)
		pls, ctg = line[0], line[1]
		pc_match = float(line[2])
		len_dict[pls] = int(line[3])
		len_dict[ctg] = int(line[4])

		if pls not in set(pls_ids):
			pls_ids.append(pls)
			pls_dict[pls] = []
		if pc_match >= 0.95:
			pls_dict[pls].append(ctg)

	predictions = {}
	#Open output file
	with open(args.detailed_output, "w") as out:
		out.write("##### general information #####\n")
		#out.write(">number of reference chromosomes: %i\n" % len(chr_dict))
		#for chr in chr_dict:
		#	out.write(chr+"\n")
		out.write(">number of reference plasmids: %i\n" % len(pls_dict))
		for pls in pls_dict:
			out.write(pls+"\n")
		out.write("\n")

		if args.plasbin != "":
			#Call function to get pb info
			print("Hello-pb")
			
			ctg_file = '../../2021-06-14__iterative_MILP_edge_novelty/output/sample_'+str(args.sample_id)+'/1.1.1/nplasmids_1/MILP/predicted_plasmids.tsv'
			string_list = read_file(ctg_file)
			for line in string_list:
				ctg = line[1:].split('\t')[1]
				pls = line.split('\t')[0]
				#print(ctg, pls)
				if pls not in predictions:
					predictions[pls] = []
				if ctg not in set(predictions[pls]):
					predictions[pls].append(ctg)
			
			eval_results(predictions, pls_dict, len_dict, len_th, out)

		if args.plasbin_flow == 1:
			#Call function to get pbf info
			
			#ctg_file = out+'/predicted_plasmids.tsv'
			ctg_file = out+'/pbf_bins.out'
			string_list = read_file(ctg_file)
			for line in string_list:
				#ctg = line[1:].split('\t')[1]
				#pls = line.split('\t')[0]
				#print(ctg, pls)
				#if pls not in predictions:
				#	predictions[pls] = []
				#if ctg not in set(predictions[pls]):
				#	predictions[pls].append(ctg)

				pls = line.split('\t')[0]
				predictions[pls] = []

				ctg_list = line.split('\t')[-1]
				ctg_list = ctg_list.split(',')
				for ctg in ctg_list:
					ctg = ctg.split(':')[0]
					predictions[pls].append(ctg)
			
			eval_results(predictions, pls_dict, len_dict, len_th, out)

		elif args.greedy != "":
			#Call function to get hy info
			print("Hello-hy")
			predictions = {}
			ctg_file = '../../../data/HyAsP_old/sample_'+str(args.sample_id)+'/putative_plasmid_contigs.fasta'
			string_list = read_file(ctg_file)
			for line in string_list:
				if line[0] == '>':
					ctg = line[1:].split('|')[0]
					pls = line.split('_')[-1]
					#print(ctg, pls)
					if pls not in predictions:
						predictions[pls] = []
					if ctg not in set(predictions[pls]):
						predictions[pls].append(ctg)
			
			eval_results(predictions, pls_dict, len_dict, len_th, out)

		elif args.gplas != "":
			#Call function to get gp info
			print("Hello-gp")
			predictions = {}
			ctg_file = '../../../data/gplas_res/sample_'+str(args.sample_id)+'_bins.tab'
			string_list = read_file(ctg_file)
			for line in string_list[1:]:
				line = line.split(' ')
				ctg = line[0]
				pls = line[1]
				#print(ctg, pls)
				if pls not in predictions:
					predictions[pls] = []
				if ctg not in set(predictions[pls]):
					predictions[pls].append(ctg)
			
			eval_results(predictions, pls_dict, len_dict, len_th, out)

		elif args.mob_recon != "":
			#Call function to get mob_recon info
			print("Hello-mob")
			predictions = {}
			ctg_file = '../../../data/mob_recon_res/sample_'+str(args.sample_id)+'/contig_report.txt'
			string_list = read_file(ctg_file)
			for line in string_list[1:]:
				line = line.split('\t')
				pls = line[1]
				ctg = line[2].split('|')[1].split('_')[0]
				length = line[3]
				#print(ctg, pls)
				if pls != 'chromosome':
					if pls not in predictions:
						predictions[pls] = []
					if ctg not in set(predictions[pls]):
						predictions[pls].append(ctg)
			
			eval_results(predictions, pls_dict, len_dict, len_th, out)
		
		elif args.plasmidspades != "":
			print("Hello-psp")
			predictions = {}
			ctg_file = args.ground_truth + '/spades_ctg_bin.tsv'
			string_list = read_file(ctg_file)
			for line in string_list[1:]:
				line = line.split('\t')
				#print(line)
				pls, ctg = line[0], line[1]
				pc_match = float(line[2])
				len_dict[pls] = int(line[3])
				len_dict[ctg] = int(line[4])

				if pls not in predictions:
					predictions[pls] = []
				if ctg not in set(predictions[pls]):
					predictions[pls].append(ctg)
			
			eval_results(predictions, pls_dict, len_dict, len_th, out)
			


main()			


