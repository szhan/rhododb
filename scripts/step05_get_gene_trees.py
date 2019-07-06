from os import listdir
from collections import OrderedDict
from random import randint
import re


#seq_type = "nt"
seq_type = "aa"

home_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/"
in_dir = home_dir + "gene_msas/"
in_suffix = "." + seq_type + ".phy"
out_dir = home_dir + "gene_trees/" + seq_type + "/"

raxml_partition_file = home_dir + "species_trees/aa_reduced_gamma/part.raxml.aa_reduced_gamma.txt"
partitionfinder_config_file = home_dir + "partitions/aa_reduced_gamma/partition_finder.cfg"

raxml_exe = "/data/home/szhan/projects/misc/1kp/bin/standard-RAxML-8.2.11/raxmlHPC-PTHREADS-AVX"
nbr_cpus = 2
nbr_bootstraps = 100


def parse_partitionfinder_config_file(in_file):
	part_dict = OrderedDict()	# gene => partition
	with open(in_file, 'rU') as in_fh:
		for line in in_fh:
			# e.g., psaC = 1-112;
			line = line.strip().replace(' ', '')
			if re.search(r'\=\d+\-\d+', line):
				(gene, part) = line.split("=")
				part = part.replace(';', '')
				part_dict[gene] = part
	return part_dict


def parse_raxml_partition_file(in_file):
	model_dict = OrderedDict()	# partition => model
	with open(in_file, 'rU') as in_fh:
		for line in in_fh:
			# e.g., FLU, Subset1 = 1-112, 24038-24525
			line = line.strip().replace(' ', '')
			if "=" in line:
				(left, right) = line.split("=")
				model = left.split(",")[0]
				parts = right.split(",")
				for part in parts:
					model_dict[part] = model
	return model_dict


def generate_random_seed():
	return "".join([str(randint(1,9)) for _ in range(5)])


genes = [f.replace(in_suffix, '') for f in listdir(in_dir) if f.endswith(in_suffix)]

if seq_type == "nt":
	raxml_models = ['GTRGAMMA', 'GTRGAMMAI']	# nt
elif seq_type == "aa":
	raxml_models = ['PROTGAMMA']	# aa
	gene_part_dict = parse_partitionfinder_config_file(partitionfinder_config_file)
	part_model_dict = parse_raxml_partition_file(raxml_partition_file)
	gene_model_dict = OrderedDict()
	for gene, part in gene_part_dict.items():
		gene_model_dict[gene] = part_model_dict[part]

for gene in genes:
	aln_file = in_dir + gene + in_suffix
	for model in raxml_models:
		out_prefix = gene + "_" + model
		if seq_type == "aa":
			model += gene_model_dict[gene]
		p_seed = generate_random_seed()
		x_seed = generate_random_seed()
		raxml_cmd = " ".join([raxml_exe, "-T", str(nbr_cpus), "-f a",\
					"-p", p_seed, "-x", x_seed, "-#", str(nbr_bootstraps),\
					"-s", aln_file, "-w", out_dir, "-n", out_prefix, "-m", model])
		print raxml_cmd

