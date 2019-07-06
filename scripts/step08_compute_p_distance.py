from os import listdir
from collections import OrderedDict

import numpy as np
from Bio import SeqIO


def get_p_distance(s1, s2):
	nbr_sites = 0
	nbr_mismatches = 0
	for c1, c2 in zip(s1, s2):
		if c1 == '-' or c2 == '-':
			pass
		else:
			nbr_sites += 1
			if c1 != c2:
				nbr_mismatches += 1
	p_dist = float(nbr_mismatches) / float(nbr_sites)
	return p_dist


def compute_p_distance(in_file):
	p_dist_list = list()
	a_len_list = list()
	seqs = OrderedDict()
	with open(in_file, 'r') as in_fh:
		for record in SeqIO.parse(in_fh, 'fasta'):
			seq_str = str(record.seq)
			seq_len = len(seq_str)
			seqs[record.id] = seq_str
			a_len_list.append(seq_len)
	id_seqs = seqs.keys()
	nbr_seqs = len(seqs)
	for i in range(nbr_seqs - 1):
		for j in range(i + 1, nbr_seqs):
			if id_seqs[i] != id_seqs[j]:
				seq_1 = seqs[id_seqs[i]]
				seq_2 = seqs[id_seqs[j]]
				p_dist = get_p_distance(seq_1, seq_2)
				p_dist_list.append(p_dist)
	median_p_dist = np.median(p_dist_list)
	assert np.mean(a_len_list) == a_len_list[0],\
		"ERROR: Alignment length not consistent in {}".format(in_file)
	return median_p_dist, a_len_list[0]


in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/"
in_msa_dir = in_base_dir + "gene_msas/"
out_result_file = in_base_dir + "p_distance/p_dist.txt"

in_suffix_aa = ".aa.aln"
in_suffix_nt = ".nt.aln"
genes = [x.replace(in_suffix_aa, '') for x in listdir(in_msa_dir) if x.endswith(in_suffix_aa)]

with open(out_result_file, 'w') as out_fh:
	out_fh.write(",".join(['gene',\
				'median_p_dist_aa', 'a_len_aa',\
				'median_p_dist_nt', 'a_len_nt']) + "\n")
	for gene in genes:
		in_msa_prefix = in_msa_dir + gene
		in_aa_file = in_msa_prefix + in_suffix_aa
		in_nt_file = in_msa_prefix + in_suffix_nt
		(median_p_dist_aa, a_len_aa) = compute_p_distance(in_aa_file)
		(median_p_dist_nt, a_len_nt) = compute_p_distance(in_nt_file)
		out_fh.write(",".join([gene,\
				str(median_p_dist_aa), str(a_len_aa),\
				str(median_p_dist_nt), str(a_len_nt)]) + "\n")


