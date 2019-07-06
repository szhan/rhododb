from os import listdir, system
from Bio import SeqIO
import dendropy


def clean_taxon_label(taxon_label):
	return taxon_label.replace('.', '_').replace('-', '_')


def prepare_mrbayes_run(gene_name, msa_file, gene_tree_file, species_tree_file, reroot_taxa, base_dir,\
			path_mrbayes_exe, path_mrbayes_script_template, tree_schema="newick"):
	out_mrbayes_dir = base_dir + "/" + gene_name + "/"
	out_nxs_file = out_mrbayes_dir + gene_name + ".clean.nxs"
	out_tree_file = out_mrbayes_dir + gene_name + ".pruned.tre"
	out_script_file = out_mrbayes_dir + gene_name + ".mrbayes.script"
	
	# prune taxa found in species tree but not gene tree
	tree_species = dendropy.Tree.get(path=species_tree_file, preserve_underscores=True, schema=tree_schema)
	
	# set branch lengths to null so as to obtain topology only
	# also, clean taxon labels while at it
	for node in tree_species.postorder_node_iter():
		node.label = None
		node.edge_length = None
		if node.is_leaf():
			node.taxon.label = clean_taxon_label(node.taxon.label)
	
	# reroot species tree and prune out the outgroup taxa
	reroot_mrca = tree_species.mrca(taxon_labels=reroot_taxa)
	tree_species.reroot_at_node(reroot_mrca, update_bipartitions=False)
	tree_species.prune_taxa_with_labels(reroot_taxa)
	assert tree_species.is_rooted == True, 'ERROR: processed species tree is not rooted!'
	taxa_species = set([node.taxon.label for node in tree_species.leaf_nodes()])
	
	tree_gene = dendropy.Tree.get(path=gene_tree_file, preserve_underscores=True, schema=tree_schema)
	taxa_gene = set([clean_taxon_label(node.taxon.label) for node in tree_gene.leaf_nodes()])
	
	# prune extra taxa from species tree
	taxa_to_prune = list(taxa_species - taxa_gene)
	if len(taxa_to_prune) > 0:
		tree_species.prune_taxa_with_labels(taxa_to_prune)
	
	# get updated set of taxa from pruned species tree
	taxa_species = set([node.taxon.label for node in tree_species.leaf_nodes()])
	
	# output files needed for mrbayes run
	system("mkdir " + out_mrbayes_dir)
	
	# prepare NEXUS file
	with open(msa_file, 'r') as in_fh, open(out_nxs_file, 'w') as out_fh:
		seqs = dict()
		for record in SeqIO.parse(in_fh, 'fasta'):
			clean_record_id = clean_taxon_label(record.id)
			if clean_record_id in taxa_species:
				seqs[clean_record_id] = str(record.seq)
		
		nbr_taxa = len(taxa_species)
		nbr_char = len(seqs.values()[0])
		
		out_fh.write("begin data;\n")
		out_fh.write(" dimensions ntax=" + str(nbr_taxa) + " " + "nchar=" + str(nbr_char) + ";\n")
		out_fh.write(" format datatype=dna gap=-;\n")
		out_fh.write("  matrix\n")
		for id, seq in seqs.items():
			out_fh.write(id + "\t" + seq + "\n")
		out_fh.write("	;\n")
		out_fh.write("end;\n")
	
	# prepare run file
	tree_newick = tree_species.as_string(schema=tree_schema).replace("'", '').strip()
	with open(path_mrbayes_script_template, 'r') as in_fh, open(out_script_file, 'w') as out_fh:
		for line in in_fh:
			line = line.strip()
			if line.endswith('TREE_NEWICK'):
				line = line.replace('TREE_NEWICK', tree_newick)
			elif line.endswith('NXS_FILE;'):
				line = line.replace('NXS_FILE', out_nxs_file)
			elif line.endswith('TREE_FILE;'):
				line = line.replace('TREE_FILE', out_tree_file)
			out_fh.write(line + "\n")
	
	mrbayes_cmd = " ".join(["cd", out_mrbayes_dir, ";", path_mrbayes_exe, out_script_file])
	
	return mrbayes_cmd


# prepare traitrates files, gene by gene
in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/"
in_suffix = ".nt.aln.trimmed"
in_species_tree_file = "/data/home/szhan/projects/misc/cp-red/analysis/species_trees_concat/RAxML_bestTree.concat_cds_aa"
out_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/mb_trees_dups/"

mrbayes_exe = "/data/home/szhan/bin/mrbayes-3.2.6/src/mb"
mrbayes_script_template = "/data/home/szhan/projects/misc/cp-red/scripts/mrbayes.script.template"

in_msa_dir = in_base_dir + "gene_msas/"
in_tree_dir = in_base_dir + "gene_trees/"

reroot_taxa = [	'Cyanidioschyzon_merolae_Cyanidiales_Cyanidiaceae_NC_004799_1',
		'Cyanidium_sp__Cyanidiales_Cyanidiaceae_KJ569775_1',
		'Cyanidium_caldarium_Cyanidiales_Cyanidiaceae_NC_001840_1',
		'Galdieria_partita_Cyanidiales_Galdieriaceae_GPAL',
		'Galdieria_sulphuraria_Cyanidiales_Galdieriaceae_NC_024665_1'
		]

#gene_names = [x.replace(in_suffix, '') for x in listdir(in_msa_dir) if x.endswith(in_suffix)]
gene_names = ['atpB', 'cpcA', 'cpcB', 'psaA', 'psbA', 'rps3']

for in_gene_name in gene_names:
	in_msa_file = in_msa_dir + in_gene_name + in_suffix
	in_gene_tree_file = in_tree_dir + "RAxML_bestTree." + in_gene_name + "_nt"
	mrbayes_cmd = prepare_mrbayes_run(gene_name=in_gene_name, msa_file=in_msa_file,
						gene_tree_file=in_gene_tree_file, species_tree_file=in_species_tree_file,
						reroot_taxa=reroot_taxa, base_dir=out_base_dir,
						path_mrbayes_exe=mrbayes_exe, path_mrbayes_script_template=mrbayes_script_template)
	print mrbayes_cmd

