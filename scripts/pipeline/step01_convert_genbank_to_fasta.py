import sys

from Bio import SeqIO
from Bio.Seq import Seq


'''
Extracts CDS features and their sequences from GenBank entries,\
prints them out in a FastA file with description lines formatted.
'''

home_dir = ""
data_dir = home_dir + "analysis2/data/"				# CHANGE
info_file = home_dir + "metadata/accessions_analysis2.txt"	# CHANGE


gb_file = data_dir + "sequence.gb"
nt_file = data_dir + "sequence.cds.fa"

info = dict()
for line in open(info_file, 'rU'):
	temp = line.strip().split("\t")
	organism = temp[0]
	order = temp[4]
	family = temp[5]
	accession = temp[6]
	info[accession] = [organism, order, family]


with open(gb_file, 'rU') as gb_fh,\
	open(nt_file, 'w') as nt_fh:
	for record in SeqIO.parse(gb_fh, "genbank"):
		if record.id not in info:
			continue
		(organism, order, family) = info[record.id]
		for feature in record.features:
			if feature.type == "CDS":
				if 'gene' in feature.qualifiers:
					gene = feature.qualifiers['gene'][0].replace('-', '_')
					protein_id = feature.qualifiers['protein_id'][0]
					sequence = str(feature.extract(record.seq)).strip('[]')
					desc_line = "-".join([organism, order, family, gene, protein_id, record.id])
					nt_fh.write(">" + desc_line + "\n")
					nt_fh.write(sequence + "\n")


