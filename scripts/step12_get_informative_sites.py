#/data/home/szhan/bin/anaconda3/bin/python
from os import listdir
from amas import AMAS


in_msa_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/gene_msas/"
in_suffix = ".nt.aln"
in_data_type = "dna"
out_stats_file = in_data_type + ".stats.txt"

genes = [f.replace(in_suffix, "") for f in listdir(in_msa_dir) if f.endswith(in_suffix)]

with open(out_stats_file, "w") as fh:
	for gene in genes:
		msa_file = in_msa_dir + gene + in_suffix
		msa = AMAS.MetaAlignment(in_files=[msa_file], data_type=in_data_type, in_format="fasta", cores=1)
		stats = msa.get_summaries()
		stats[1][0][0] = stats[1][0][0].replace(in_suffix, "")
		fh.write(",".join(stats[1][0]) + "\n")


