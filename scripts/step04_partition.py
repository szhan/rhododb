from os import listdir
from collections import OrderedDict
from Bio import SeqIO, AlignIO


home_dir = "/data/home/szhan/projects/misc/cp-red/"
aln_dir = home_dir + "analysis2/gene_msas/"		# CHANGE
part_dir = home_dir + "analysis2/partitions/"		# CHANGE

in_aa_suffix = '.aa.aln'
in_nt_suffix = '.nt.aln'
out_aa_suffix = '.aa.phy'
out_nt_suffix = '.nt.phy'


genes = [x.replace(in_aa_suffix, '') for x in listdir(aln_dir) if x.endswith(in_aa_suffix)]
for gene in genes:
	in_aa_file = aln_dir + gene + in_aa_suffix
	in_nt_file = aln_dir + gene + in_nt_suffix
	out_aa_file = aln_dir + gene + out_aa_suffix
	out_nt_file = aln_dir + gene + out_nt_suffix
	# convert to phylip
	SeqIO.convert(in_aa_file, "fasta", out_aa_file, "phylip-relaxed")
	SeqIO.convert(in_nt_file, "fasta", out_nt_file, "phylip-relaxed")


def concatenate_msa(genes, in_suffix, out_concat_file, out_part_file):
	concat_dict = OrderedDict()
	concat_part = OrderedDict()
	
	for gene in genes:
		in_file = aln_dir + gene + in_suffix
		msa_aln = AlignIO.read(open(in_file, 'rU'), 'fasta')
		for record in msa_aln:
			concat_dict[record.id] = ''
	
	total_len = 0
	for gene in genes:
		in_file = aln_dir + gene + in_suffix
		msa_aln = AlignIO.read(open(in_file, 'rU'), 'fasta')
		msa_len = msa_aln.get_alignment_length()
		# first, add sequences in MSA
		visited = list()
		for record in msa_aln:
			concat_dict[record.id] += str(record.seq)
			visited.append(record.id)
		# second, add hyphens in sequences not in MSA
		for id in concat_dict.keys():
			if id not in visited:
				concat_dict[id] += '-' * msa_len
		# get partitioning scheme
		concat_part[gene] = str(total_len + 1) + '-' + str(total_len + msa_len)
		total_len += msa_len
	
	with open(out_concat_file, 'w') as out_fh:
		out_fh.write(" " + str(len(concat_dict)) + " " + str(total_len) + "\n")
		for id, seq in concat_dict.items():
			out_fh.write(id + "\t" + seq + "\n")
	
	with open(out_part_file, 'w') as out_fh:
		for id in concat_part.keys():
			out_fh.write(id + " = " + concat_part[id] + ";\n")


out_aa_concat_file = part_dir + 'concat.aa.phy'
out_nt_concat_file = part_dir + 'concat.nt.phy'
out_aa_part_file = part_dir + 'part.aa.txt'
out_nt_part_file = part_dir + 'part.nt.txt'
concatenate_msa(genes, in_aa_suffix, out_aa_concat_file, out_aa_part_file)
concatenate_msa(genes, in_nt_suffix, out_nt_concat_file, out_nt_part_file)


