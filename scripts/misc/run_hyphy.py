from os import listdir, system
import dendropy


# Remove periods and hyphens in taxon names, as HyPhy cannot handle them
def clean_taxon_label(taxon_label):
	return taxon_label.replace('.', '_').replace('-', '_')


def prepare_hyphy_run(gene_name, msa_file, gene_tree_file, species_tree_file, base_dir,\
			path_hyphy_exe, path_hyphy_batch_file, tree_schema="newick"):
	out_hyphy_dir = base_dir + "/" + gene_name + "/"
	out_msa_file = out_hyphy_dir + gene_name + ".clean.aln"
	out_tree_file = out_hyphy_dir + gene_name + ".pruned.tre"
	out_result_file = out_hyphy_dir + gene_name + ".result.txt"
	
	# get taxa found in species tree but not gene tree
	tree_gene = dendropy.Tree.get(path=gene_tree_file, preserve_underscores=True, schema=tree_schema)
	taxa_gene = set([node.taxon.label for node in tree_gene.leaf_nodes()])
	
	tree_species = dendropy.Tree.get(path=species_tree_file, preserve_underscores=True, schema=tree_schema)
	taxa_species = set([node.taxon.label for node in tree_species.leaf_nodes()])
	
	assert len(taxa_species) >= len(taxa_gene)
	taxa_to_prune = list(taxa_species - taxa_gene)
	
	# prune extra taxa from species tree
	if len(taxa_to_prune) > 0:
		tree_species.prune_taxa_with_labels(taxa_to_prune)
	
	# set branch lengths to null so as to obtain topology only
	for node in tree_species.postorder_node_iter():
		node.label = None
		node.edge_length = None
		if node.is_leaf():
			node.taxon.label = clean_taxon_label(node.taxon.label)
	
	# replace quote characters around taxon labels
	tree_str = tree_species.as_string(schema=tree_schema).replace("'", '').strip()
	
	# output files needed for hyphy run
	system("mkdir " + out_hyphy_dir)
	
	with open(out_tree_file, 'w') as out_fh:
		out_fh.write(tree_str)
	
	with open(msa_file, 'r') as in_fh,\
		open(out_msa_file, 'w') as out_fh:
			for line in in_fh:
				line = line.strip()
				if line.startswith('>'):
					line = clean_taxon_label(line[1:])
					out_fh.write(">" + line + "\n")
				else:
					out_fh.write(line + "\n")
	
	hyphy_cmd = " ".join(["cd", out_hyphy_dir, ";"])
	hyphy_cmd += " ".join([path_hyphy_exe, path_hyphy_batch_file, "Universal", out_msa_file, out_tree_file, "All", ">", out_result_file])
	
	return hyphy_cmd


# prepare hyphy files, gene by gene
in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/"
in_suffix = ".nt.aln.trimmed"
in_species_tree_file = "/data/home/szhan/projects/misc/cp-red/analysis/species_trees_concat/RAxML_bestTree.concat_cds_aa"
out_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/hyphy_dups/"

num_cpus = 16
mpirun_exe = "mpirun"
hyphy_exe = "/data/home/szhan/projects/misc/cp-red/bin/hyphy-2.3.12_install/bin/HYPHYMPI"
hyphympi_exe = " ".join([mpirun_exe, "-np", str(num_cpus), hyphy_exe])
hyphy_batch_file = "/data/home/szhan/projects/misc/cp-red/bin/hyphy-2.3.12_install/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf"

#gene_names = [x.replace(in_suffix, '') for x in listdir(in_base_dir + "gene_msas/") if x.endswith(in_suffix)]
gene_names = ['atpB', 'cpcA', 'cpcB', 'psaA', 'psbA', 'rps3']

for in_gene_name in gene_names:
	in_msa_file = in_base_dir + "gene_msas/" + in_gene_name + in_suffix
	in_gene_tree_file = in_base_dir + "gene_trees/" + "RAxML_bestTree." + in_gene_name + "_nt"
	hyphy_cmd = prepare_hyphy_run(gene_name=in_gene_name, msa_file=in_msa_file,\
					gene_tree_file=in_gene_tree_file, species_tree_file=in_species_tree_file, base_dir=out_base_dir,\
					path_hyphy_exe=hyphympi_exe, path_hyphy_batch_file=hyphy_batch_file)
	print hyphy_cmd


