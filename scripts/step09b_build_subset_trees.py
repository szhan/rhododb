from collections import OrderedDict
from itertools import combinations
from os import listdir
import sys

from Bio import SeqIO
from numpy import mean


seq_type = sys.argv[1]
subset_size = int(sys.argv[2])

assert seq_type in ['aa', 'nt']

in_base_dir = ""
in_file = in_base_dir + "metadata/accessions_analysis2.txt"
in_dir = in_base_dir + "analysis2/gene_msas/"
out_dir = in_base_dir + "analysis2/subsets/" + seq_type + "/"

in_suffix = "." + seq_type + ".aln"

bin_dir = "/data/home/szhan/projects/misc/1kp/bin/"
raxml_exe = bin_dir + "standard-RAxML-8.2.11/raxmlHPC-PTHREADS-AVX"


model = "PROTGAMMAWAG" if seq_type == "aa" else "GTRGAMMA"

genes = ['apcE', 'gltB', 'rpoc1', 'rpoc2', 'rpoB']
subsets = combinations(genes, subset_size)

for subset in subsets:
	out_prefix = "_".join(subset) + ".subset." + seq_type
	out_aln_file = out_dir + "/" + out_prefix + ".aln"
	out_part_file = out_dir + "/" + out_prefix + ".part"

	all_taxa = OrderedDict()

	with open(in_file, 'r') as fh:
		for line in fh:
			info = line.split("\t")
			species = info[0]
			order = info[4]
			family = info[5]
			accession = info[6]
			# e.g., Kuetzingia_canaliculata-Ceramiales-Rhodomelaceae-MF101449.1
			member = "-".join([species, order, family, accession])
			all_taxa[member] = ""

	start = 1
	parts = list()

	for gene in subset:
		in_taxa = set(all_taxa.keys())
		in_gene_file = in_dir + gene + in_suffix

		with open(in_gene_file, 'r') as fh:
			# concatenate sequence
			aln_len_list = list()
			for record in SeqIO.parse(fh, 'fasta'):
				seq_id = record.id
				seq_seq = str(record.seq)
				aln_len_list.append(len(seq_seq))
				all_taxa[seq_id] += seq_seq
				if seq_id in in_taxa:
					in_taxa.remove(seq_id)

			# check alignment length
			aln_len = aln_len_list[0]
			assert mean(aln_len_list) == aln_len,\
				"ERROR: Alignment length is inconsistent in {}".format(in_gene_file)

			# add gap character to missing entries
			for id in in_taxa:
				all_taxa[id] += '-' * aln_len

			# print partitions
			bounds = str(start) + "-" + str(start + aln_len - 1)
			sub_model = 'CPREV' if seq_type == 'aa' else 'DNA'
			part_info = ' '.join([sub_model + ',', gene, '=', bounds])
			parts.append(part_info)
			start += aln_len

	with open(out_part_file, 'w') as fh:
		for part in parts:
			fh.write(part + "\n")

	with open(out_aln_file, 'w') as fh:
		for id in all_taxa.keys():
			seq_len_wo_gaps = len(all_taxa[id].replace('-', ''))
			if seq_len_wo_gaps > 0:
				fh.write(">" + id + "\n")
				fh.write(all_taxa[id] + "\n")

	raxml_pars = ["-T 8", "-f a", "-p 56241"]#, "-x 15873", "-# 100"]
	raxml_cmd = " ".join(	[raxml_exe] + raxml_pars + ["-m", model] +\
				["-w", out_dir, "-n" , out_prefix,\
				"-s", out_aln_file, "-q", out_part_file])
	print raxml_cmd
