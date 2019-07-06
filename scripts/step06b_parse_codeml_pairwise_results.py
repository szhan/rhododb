import re
from os import listdir
from collections import OrderedDict

import numpy as np


def parse_codeml_pairwise_result_file(in_file):
	nbr_entries = None
	pairwise_matrix = OrderedDict()
	with open(in_file, 'r') as in_fh:
		for line in in_fh:
			line = line.rstrip()
			match = re.match(r'\s+(\d+)', line)
			if match:
				nbr_entries = int(match.group(1))
			else:
				# keep elements from half a symmetric matrix
				elements = re.split(r"\s+", line)
				if len(elements) == 1:
					pairwise_matrix[elements[0]] = []
				else:
					pairwise_matrix[elements[0]] = elements[1:]
	assert nbr_entries == len(pairwise_matrix),\
		"ERROR: pairwise matrix contains {} rows, when expecting {}".format(nbr_entries, len(pairwise_matrix))
	return pairwise_matrix


def compute_dNdS_statistics(in_file_dN, in_file_dS):
	codeml_dN = parse_codeml_pairwise_result_file(in_file_dN)
	codeml_dS = parse_codeml_pairwise_result_file(in_file_dS)

	# check that two matrices correspond
	assert len(codeml_dN) == len(codeml_dS), "ERROR: dN and dS matices got different number of rows, {} and {}".format(len(codeml_dN), len(codeml_dS))
	assert codeml_dN.keys() == codeml_dS.keys(), "ERROR: dN and dS matrices got different rows"

	all_dN = np.array([])
	all_dS = np.array([])
	all_dNdS = np.array([])

	nbr_rows = len(codeml_dN)
	nbr_elements = (nbr_rows * (nbr_rows - 1)) / 2
	row_names = codeml_dN.keys()

	for i in range(len(codeml_dN)):
		row_dN = np.array(codeml_dN[row_names[i]], dtype='float')
		row_dS = np.array(codeml_dS[row_names[i]], dtype='float')
		assert len(row_dN) == len(row_dS),\
			"ERROR: number of elements in row {} not identical in dN and dS".format(str(i))

		if len(row_dN) == 0:
			continue

		index_valid_dS = [i for i in range(len(row_dS)) if row_dS[i] > 0 and not np.isnan(row_dS[i])]
		row_dN = np.array([row_dN[i] for i in index_valid_dS])
		row_dS = np.array([row_dS[i] for i in index_valid_dS])
		row_dNdS = np.divide(row_dN, row_dS)

		all_dN = np.append(all_dN, row_dN)
		all_dS = np.append(all_dS, row_dS)
		all_dNdS = np.append(all_dNdS, row_dNdS)

	assert len(all_dN) == len(all_dS) == len(all_dNdS),\
		"ERROR: number of elements in dN, dS, and dNdS arrays is different"

	first_quartile_dS = np.percentile(all_dS, 25)
	median_dS = np.median(all_dS)
	third_quartile_dS = np.percentile(all_dS, 75)

	first_quartile_dN = np.percentile(all_dN, 25)
	median_dN = np.median(all_dN)
	third_quartile_dN = np.percentile(all_dN, 75)

	first_quartile_dNdS = np.percentile(all_dNdS, 25)
	median_dNdS = np.median(all_dNdS)
	third_quartile_dNdS = np.percentile(all_dNdS, 75)

	dNdS_stats = [	nbr_rows, nbr_elements, len(all_dN),
			first_quartile_dS, median_dS, third_quartile_dS,
			first_quartile_dN, median_dN, third_quartile_dN,
			first_quartile_dNdS, median_dNdS, third_quartile_dNdS
			]

	return dNdS_stats


in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/"
in_suffix = ".nt.aln"
in_result_dir = in_base_dir + "/codeml_pairwise/"
out_result_file = in_result_dir + "/codeml_pairwise_results.txt"

#method = "ML"
method = "NG"

genes = [x.replace(in_suffix, '') for x in listdir(in_base_dir + "gene_msas/") if x.endswith(in_suffix)]

with open(out_result_file, 'w') as out_fh:
	out_fh.write(",".join(['gene', 'nbr_taxa', 'nbr_comps', 'nbr_valid_comps',\
				'dS_q25', 'dS_q50', 'dS_q75',\
				'dN_q25', 'dN_q50', 'dN_q75',\
				'dNdS_q25', 'dNdS_q50', 'dNdS_q75']) + "\n")

	for gene in genes:
		gene_file_prefix = in_result_dir + "/" + gene + "/2" + method
		gene_file_dN = gene_file_prefix + ".dN"
		gene_file_dS = gene_file_prefix + ".dS"
		gene_dNdS_stats = compute_dNdS_statistics(gene_file_dN, gene_file_dS)
		out_fh.write(",".join([gene] + [str(i) for i in gene_dNdS_stats]) + "\n")


