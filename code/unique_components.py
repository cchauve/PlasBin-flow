from sys import argv
import os
import argparse
import preprocessing

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--ag", help="Path to assembly graph file")
	parser.add_argument("--comps", help="Path to components file")	
	#parser.add_argument("--gc", help="Path to GC probabilities file")
	#parser.add_argument("--c", help="Path to contig classification file")
	#parser.add_argument("--seeds", help="Path to seed contigs file")
	parser.add_argument("--out", help="Path to output dir")
	#parser.add_argument("--alpha1", nargs='?', const = 1, type=int, default = 1, help="Weight of flow term")
	#parser.add_argument("--alpha2", nargs='?', const = 1, type=int, default = 1, help="Weight of GC content term")
	#parser.add_argument("--alpha3", nargs='?', const = 1, type=int, default = 1, help="Weight of log probabilities term")	
	#parser.add_argument("--rmiter", nargs='?', const = 1, type=int, default = 50, help="Number of iterations to remove circular components")


	args = parser.parse_args()

	output_dir = args.out
	assembly_file = args.ag
	components_file = args.comps
	#gc_file = args.gc
	#class_file = args.c
	#seeds_file = args.seeds
	#alpha1 = args.alpha1
	#alpha2 = args.alpha2
	#alpha3 = args.alpha3
	#rmiter = args.rmiter

	string_list = read_file(components_file)

	unique_comps = []
	for line in string_list:
		line = line.split("\t")[1]
		line = line.split(",")
		#print(line)

		current_comp = set()
		for contig in line:
			if contig not in ['S','T','']:
				current_comp.add(contig)
		#print(current_comp)	

		current_unique_comps = []
		flag = 1	#current component is unique
		for comp in unique_comps:
			print(comp)
			if not comp.issubset(current_comp):
				current_unique_comps.append(comp)
			if current_comp.issubset(comp):
				flag = 0
		if flag == 1:
			current_unique_comps.append(current_comp)
		unique_comps = current_unique_comps

	#print(unique_comps)

	#Getting contig info
	contigs_dict = {}
	links_list = []
	contigs_dict, links_list = preprocessing.get_data(assembly_file, contigs_dict, links_list)

	unique_components = 'unique_components.csv'
	unique_components_file = open(os.path.join(output_dir, unique_components), "w")

	formatted_input = 'predicted_plasmids.tsv'
	formatted_input_file = open(os.path.join(output_dir, formatted_input), "w")
	formatted_input_file.write("#Plasmid\tContig\tLength\n")

	count = 0
	for comp in unique_comps:
		count += 1
		unique_components_file.write("Component_"+str(count)+";")
		for contig in comp:
			unique_components_file.write(contig+",")
			c_len = contigs_dict[contig]['Length']
			formatted_input_file.write("P"+str(count)+"\t"+contig+"\t"+str(c_len)+"\n")
		unique_components_file.write("\n")	
		formatted_input_file.write("\n")









			


		


