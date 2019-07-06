import re
import sys
from os import system
from collections import Counter
from Bio import SeqIO


# ran TransDecoder as follows:
#TD_DIR="/data/home/szhan/projects/misc/1kp-algae/bin/TransDecoder-3.0.0/"
#$TD_DIR"/TransDecoder.LongOrfs" -m 50 -t sequence.cds.fa
#$TD_DIR"/TransDecoder.Predict" --single_best_orf -t sequence.cds.fa


home_dir = "/data/home/szhan/projects/misc/cp-red/"
info_file = home_dir + "metadata/accessions_analysis2.txt"	# CHANGE

in_dir = home_dir + "analysis2/data/"				# CHANGE
in_nt_file = in_dir + "sequence.cds.fa.transdecoder.cds"
in_aa_file = in_dir + "sequence.cds.fa.transdecoder.pep"
out_aln_dir = home_dir + "analysis2/gene_msas/"			# CHANGE
out_tre_dir = home_dir + "analysis2/gene_trees/"		# CHANGE

system("mkdir " + out_aln_dir)
system("mkdir " + out_tre_dir)

# check manually
correct_names = {'chlI_2':'chlI', 'APC_beta':'apcB'}
min_taxa_per_cluster = 96					# CHANGE


def parse_accession_file(in_file):
	"""
	Get accession to species map.
	e.g., Galdieria_sulphuraria,Rhodophyta,Cyanidiophyceae,-,Cyanidiales,Galdieriaceae,NC_024665.1
		turned to NC_024665.1 -> Galdieria_sulphuraria
	"""
	accession_species_map = dict()
	with open(in_file, 'rU') as fh:
		for line in fh:
			info = line.split("\t")
			species = info[0]
			accession = info[6]
			accession_species_map[accession] = species
	return accession_species_map


def get_original_entry_name(desc_line):
	"""
	Parse description line
	e.g., >Gene.22020::Acrosorium_ciliolatum-Ceramiales-Rhodomelaceae-ConsOrf1-ARW59975.1-MF101411.1::g.22020::m.22020 Gene.22020::Acrosorium_ciliolatum-Ceramiales-Rhodomelaceae-ConsOrf1-ARW59975.1-MF101411.1::g.22020  ORF type:complete len:61 (+) Acrosorium_ciliolatum-Ceramiales-Rhodomelaceae-ConsOrf1-ARW59975.1-MF101411.1:1-183(+)
	"""
	return re.split(r'\s+', desc_line)[0].split('::')[1]


def reformat_gene_name(gene_name):
	""" Clean gene name. """
	gene_name = gene_name.lower()
	rd = re.compile(r"\d+$")	# e.g., rpl17
	if len(gene_name) == 4 and not rd.match(gene_name):
		gene_name = gene_name[:-1] + gene_name[-1].upper()
	return gene_name


# load sequences
# assuming that original description lines in NT and AA files should be identical
seqs_nt = dict()
with open(in_nt_file, 'rU') as nt_fh:
	for record in SeqIO.parse(nt_fh, "fasta"):
		id = get_original_entry_name(record.id)
		seqs_nt[id] = str(record.seq)

seqs_aa = dict()
with open(in_aa_file, 'rU') as aa_fh:
	for record in SeqIO.parse(aa_fh, "fasta"):
		id = get_original_entry_name(record.id)
		seqs_aa[id] = str(record.seq)
		if seqs_aa[id].endswith('*'):
			seqs_aa[id] = seqs_aa[id][:-1]


# cluster genes by gene name in description line
gene_families = dict()
for id in seqs_aa.keys():
	gene_name = id.split('-')[3]
	if gene_name in correct_names:
		gene_name = correct_names[gene_name]
	gene_name = reformat_gene_name(gene_name)
	#if gene_name.endswith('fragment'):
	#	continue
	if gene_name in gene_families:
		gene_families[gene_name] += [id]
	else:
		gene_families[gene_name] = [id]

# remove the genes of a species with duplicate entries in a gene family
for gene_family in gene_families.keys():
	all_accessions = [x.split('-')[5] for x in gene_families[gene_family]]
	unique_accessions = set(all_accessions)
	if len(all_accessions) != len(unique_accessions):
		print gene_family	# these got filtered out initially, but added back
		#del gene_families[gene_family]

# print each gene family into separate file
accession_species_map = parse_accession_file(info_file)
for gene_family in gene_families:
	if len(gene_families[gene_family]) >= min_taxa_per_cluster:
		out_aln_prefix = out_aln_dir + '/' + gene_family
		out_tre_prefix = out_tre_dir + '/' + gene_family
		
		out_nt_fa_file = out_aln_prefix + '.nt.fa'
		out_aa_fa_file = out_aln_prefix + '.aa.fa'
		out_nt_aln_file = out_aln_prefix + '.nt.aln'
		out_aa_aln_file = out_aln_prefix + '.aa.aln'
		#out_nt_trim_file = out_nt_aln_file + '.trimmed'
		#out_aa_trim_file = out_aa_aln_file + '.trimmed'
		
		out_nt_tre_file = out_tre_prefix + '.nt.tre'
		out_aa_tre_file = out_tre_prefix + '.aa.tre'
		
		# print to fasta file
		with open(out_nt_fa_file, 'w') as out_nt_fh,\
			open(out_aa_fa_file, 'w') as out_aa_fh:
			for member in gene_families[gene_family]:
				# remove gene name in description line
				# e.g., Hildenbrandia_rivularis-Hildenbrandiales-Hildenbrandiaceae-rpl13-AOM67222.1-KX284723.1
				info = member.split('-')
				species_name = accession_species_map[info[5]]
				order_name = info[1]
				family_name = info[2]
				accession = info[5]
				
				reformatted_desc_line = "-".join([species_name, order_name, family_name, accession])
				out_nt_fh.write(">" + reformatted_desc_line + "\n")
				out_nt_fh.write(seqs_nt[member] + "\n")
				out_aa_fh.write(">" + reformatted_desc_line + "\n")
				out_aa_fh.write(seqs_aa[member] + "\n")

