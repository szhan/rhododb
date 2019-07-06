from os import listdir
from collections import OrderedDict
import sys

import numpy as np
from Bio import SeqIO


def parse_fasta_file(in_file):
	seqs = OrderedDict()
	with open(in_file, 'r') as in_fh:
		for record in SeqIO.parse(in_fh, 'fasta'):
			seqs[record.id] = str(record.seq)
	return seqs


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
	p_dist = float(nbr_mismatches) / float(nbr_sites) if nbr_sites > 0 else 1
	return p_dist


def compute_p_distance(in_seqs, left_coord, right_coord):
	p_dist_list = list()
	id_seqs = in_seqs.keys()
	nbr_seqs = len(in_seqs)
	for i in range(nbr_seqs - 1):
		for j in range(i + 1, nbr_seqs):
			if id_seqs[i] != id_seqs[j]:
				seq_1 = in_seqs[id_seqs[i]][left_coord:right_coord]
				seq_2 = in_seqs[id_seqs[j]][left_coord:right_coord]
				p_dist = get_p_distance(seq_1, seq_2)
				p_dist_list.append(p_dist)
	median_p_dist = np.median(p_dist_list)
	mean_p_dist = np.mean(p_dist_list)
	return([median_p_dist, mean_p_dist])


in_msa_file = sys.argv[1]
in_seqs = parse_fasta_file(in_msa_file)
msa_len = len(in_seqs.values()[0])

window_size = 30
for i in range(msa_len - window_size):
	[median_p_dist, mean_p_dist] = compute_p_distance(in_seqs, i, i + window_size)
	print ",".join([str(i), str(median_p_dist), str(mean_p_dist)])


