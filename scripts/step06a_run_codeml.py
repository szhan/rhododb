from os import listdir, system
import dendropy


def clean_taxon_label(taxon_label):
	return taxon_label.replace('.', '_').replace('-', '_')


def prepare_codeml_run(gene_name, msa_file, gene_tree_file, species_tree_file, base_dir,\
			path_codeml_exe, path_codeml_ctrl, tree_schema="newick"):
	out_codeml_dir = base_dir + "/" + gene_name + "/"
	out_msa_file = out_codeml_dir + gene_name + ".aln"
	out_tree_file = out_codeml_dir + gene_name + ".pruned.tre"
	out_control_file = out_codeml_dir + "codeml.ctl"
	out_result_file = out_codeml_dir + gene_name + ".result.txt"
	
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
	
	# output files needed for codeml run
	system("mkdir " + out_codeml_dir)
	
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
	
	with open(path_codeml_ctrl, 'r') as in_fh,\
		open(out_control_file, 'w') as out_fh:
			for line in in_fh:
				line = line.strip()
				if line.endswith('MSA_FILE'):
					line = line.replace('MSA_FILE', out_msa_file)
				elif line.endswith('TREE_FILE'):
					line = line.replace('TREE_FILE', out_tree_file)
				elif line.endswith('RESULT_FILE'):
					line = line.replace('RESULT_FILE', out_result_file)
				out_fh.write(line + "\n")
	
	codeml_cmd = " ".join(["cd", out_codeml_dir, ";", path_codeml_exe])
	
	return codeml_cmd


# prepare codeml files, gene by gene
in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis2/"
in_suffix = ".nt.aln"
in_species_tree_file = in_base_dir + "species_trees/nt_reduced/RAxML_bestTree.concat_nt_reduced"
out_base_dir = in_base_dir + "codeml_pairwise/"

codeml_exe = "/data/home/szhan/projects/misc/cp-red/bin/paml4.9h/bin/codeml"
codeml_ctrl_template = in_base_dir + "../scripts/codeml.ctl.template.pairwise"

genes = [x.replace(in_suffix, '') for x in listdir(in_base_dir + "gene_msas/") if x.endswith(in_suffix)]

for gene in genes:
	in_msa_file = in_base_dir + "gene_msas/" + gene + in_suffix
	in_gene_tree_file = in_base_dir + "gene_trees/nt/RAxML_bestTree." + gene + "_GTRGAMMA"
	codeml_cmd = prepare_codeml_run(gene_name=gene, msa_file=in_msa_file,\
					gene_tree_file=in_gene_tree_file, species_tree_file=in_species_tree_file, base_dir=out_base_dir,\
					path_codeml_exe=codeml_exe, path_codeml_ctrl=codeml_ctrl_template)
	print codeml_cmd

