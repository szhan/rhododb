from os import listdir
from collections import OrderedDict
from random import randint
import sys

import numpy as np
from Bio import SeqIO


def parse_fasta_file(in_file):
	seqs = OrderedDict()
	with open(in_file, 'rU') as fh:
		for record in SeqIO.parse(fh, 'fasta'):
			seqs[record.id] = str(record.seq)
	return seqs


def slice_msa(in_seqs, left_coord, right_coord, out_file):
	with open(out_file, 'w') as fh:
		for id, seq in in_seqs.items():
			fh.write(">" + id + "\n")
			fh.write(seq[left_coord:right_coord] + "\n")


def generate_random_seed(seed_len=5):
	return "".join([str(randint(1,9)) for _ in range(seed_len)])


bin_dir = "/data/home/szhan/projects/misc/1kp/bin/"
raxml_exe = bin_dir + "standard-RAxML-8.2.11/raxmlHPC-PTHREADS-AVX"

in_gene = sys.argv[1]
in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/"
in_msa_dir = in_base_dir + "gene_msas/"
in_msa_nt_file = in_msa_dir + in_gene + ".nt.aln"
in_msa_aa_file = in_msa_dir + in_gene + ".aa.aln"
out_msa_dir = in_base_dir + "primer_design/"


in_seqs_nt = parse_fasta_file(in_msa_nt_file)
in_seqs_aa = parse_fasta_file(in_msa_aa_file)

msa_len_nt = len(in_seqs_nt.values()[0])
msa_len_aa = len(in_seqs_aa.values()[0])

window_size = 801	# multiple of 3, NT
skip_size = 90		# multiple of 3, NT

assert msa_len_nt == 3 * msa_len_aa, "ERROR: NT alignment length is not three times that of the AA alignment."
assert msa_len_nt % 3 == 0, "ERROR: NT alignment length is not a multiple of 3."
assert window_size % 3 == 0, "ERROR: Window size is not a multiple of 3."
assert skip_size % 3 == 0, "ERROR: Skip size is not a multiple of 3."

# scan NT alignment
for i in range(0, msa_len_nt - window_size, skip_size):
	out_msa_prefix = in_gene + ".nt." + str(i)
	out_msa_file = out_msa_dir + out_msa_prefix + ".aln"
	slice_msa(in_seqs_nt, i, i + window_size, out_msa_file)
	p_seed = generate_random_seed()
	x_seed = generate_random_seed()
	raxml_cmd = " ".join([raxml_exe, "-T 4", "-p", p_seed, "-s", out_msa_file, "-w", out_msa_dir, "-n" , out_msa_prefix, "-m GTRGAMMA"])
	print raxml_cmd


# scan AA alignment
window_size /= 3
skip_size /=3
for i in range(0, msa_len_aa - window_size, skip_size):
	out_msa_prefix = in_gene + ".aa." + str(i)
	out_msa_file = out_msa_dir + out_msa_prefix + ".aln"
	slice_msa(in_seqs_aa, i, i + window_size, out_msa_file)
	p_seed = generate_random_seed()
	x_seed = generate_random_seed()
	raxml_cmd = " ".join([raxml_exe, "-T 4", "-p", p_seed, "-s", out_msa_file, "-w", out_msa_dir, "-n" , out_msa_prefix, "-m PROTGAMMACPREV"])
	print raxml_cmd


