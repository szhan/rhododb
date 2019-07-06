import re
import sys
from os import system, listdir
from collections import Counter
from Bio import SeqIO


muscle_exe = "muscle3.8.31_i86linux64"


home_dir = ""
aln_dir = home_dir + "analysis2/gene_msas/"			# CHANGE
tre_dir = home_dir + "analysis2/gene_trees/"			# CHANGE

in_aa_suffix = '.aa.fa'
in_nt_suffix = '.nt.fa'
out_aa_suffix = '.aa.aln'
out_nt_suffix = '.nt.aln'


def codon_align(aa_aln_file, nt_raw_file, nt_aln_file):
	""" Takes protein alignment and fills it with corresponding codons. """
	
	aa_aln_dict = SeqIO.to_dict(SeqIO.parse(aa_aln_file, "fasta"))
	nt_raw_dict = SeqIO.to_dict(SeqIO.parse(nt_raw_file, "fasta"))
	assert len(aa_aln_dict) == len(nt_raw_dict), "ERROR: number of entries in AA and CDS files are not equal."
	
	nt_aln_dict = dict()
	for id, seq in nt_raw_dict.items():
		assert id in aa_aln_dict, "ERROR: {} is not found.".format(id)
		aa_aln_seq = str(aa_aln_dict[id].seq)
		nt_raw_seq = str(seq.seq)
		aa_raw_len = len(aa_aln_seq.replace('-', ''))	# exclude gaps
		nt_raw_len = len(nt_raw_seq)
		
		assert nt_raw_len % 3 == 0, "ERROR: CDS length of {} is not a multiple of 3.".format(id)
		assert nt_raw_len/3 == aa_raw_len or nt_raw_len/3 == aa_raw_len + 1,\
			"ERROR: Unexpected number of codons in {}".format(id)
		
		if nt_raw_len/3 == aa_raw_len + 1:
			# remove stop codon
			last_codon = nt_raw_seq[-3:]
			assert last_codon in ['TAG', 'TGA', 'TAA'], "ERROR: Last codong is not a recognized stop codon.".format(last_codon)
			nt_raw_seq = nt_raw_seq[:-3]
			nt_raw_len = len(nt_raw_seq)	# update length
		
		# replace AAs with corresponding codons
		codons = [nt_raw_seq[i:(i+3)] for i in range(0, nt_raw_len, 3)]
		nt_aln_seq = ''
		for i in range(len(aa_aln_seq)):
			if aa_aln_seq[i] == '-':
				nt_aln_seq += '---'
			else:
				nt_aln_seq += codons.pop(0)
		nt_aln_dict[id] = nt_aln_seq
	
	# print to FastA file
	with open(nt_aln_file, 'w') as out_fh:
		for id, seq in nt_aln_dict.items():
			out_fh.write(">" + id + "\n")
			out_fh.write(seq + "\n")


genes = [x.replace(in_aa_suffix, '') for x in listdir(aln_dir) if x.endswith(in_aa_suffix)]
for gene in genes:
	in_data_dir = aln_dir + gene
	in_aa_raw_file = in_data_dir + in_aa_suffix
	in_nt_raw_file = in_data_dir + in_nt_suffix
	out_aa_aln_file = in_data_dir + out_aa_suffix
	out_nt_aln_file = in_data_dir + out_nt_suffix
	muscle_aa_cmd = " ".join([muscle_exe, "-quiet", "-in", in_aa_raw_file, "-out", out_aa_aln_file])
	#system(muscle_aa_cmd)	# UNCOMMENT
	codon_align(out_aa_aln_file, in_nt_raw_file, out_nt_aln_file)


