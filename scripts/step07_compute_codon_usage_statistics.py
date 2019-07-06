from os import listdir, system
import numpy as np


def parse_codonw_output_file(in_file):
	perc_gc = list()
	gc3 = list()
	nc = list()

	with open(in_file, 'r') as in_fh:
		for line in in_fh:
			line = line.strip()
			if line.startswith('title'):
				continue
			data = line.split("\t")
			if data[8].startswith('*') or data[9].startswith('*') or data[10].startswith('*'):
				continue
			perc_gc.append(float(data[10]))
			gc3.append(float(data[9]))
			nc.append(float(data[8]))

	median_perc_gc = np.median(perc_gc)
	median_gc3 = np.median(gc3)
	median_nc = np.median(nc)

	return([median_perc_gc, median_gc3, median_nc])


codonw_exe = "/data/home/szhan/projects/misc/cp-red/bin/codonW/codonw"
codonw_args = "-all_indices -nomenu -silent"
codonw_suffix = ".nt.out"

in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/"
in_suffix = ".nt.aln"

out_results_file = in_base_dir + "codonw/codonw_results.txt"
tmp_dir = in_base_dir + "codonw/tmp/"
system("mkdir " + tmp_dir)

genes = [x.replace(in_suffix, '') for x in listdir(in_base_dir + "gene_msas/") if x.endswith(in_suffix)]

with open(out_results_file, 'w') as out_fh:
	out_fh.write(",".join(['gene', 'median_perc_gc', 'median_gc3', 'median_nc']) + "\n")
	for gene in genes:
		in_msa_file = in_base_dir + "/gene_msas/" + gene + in_suffix
		in_codonw_file = tmp_dir + gene + in_suffix
		out_codonw_file = tmp_dir + gene + codonw_suffix
		system("cp " + in_msa_file + " " + tmp_dir)

		codonw_cmd = " ".join([codonw_exe, in_codonw_file, codonw_args])
		system(codonw_cmd)

		stats = parse_codonw_output_file(out_codonw_file)
		out_fh.write(",".join([str(x) for x in [gene] + stats]) + "\n")


