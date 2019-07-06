from os import listdir, system
from Bio import SeqIO
import dendropy


def clean_taxon_label(taxon_label):
	return taxon_label.replace('.', '_').replace('-', '_')


def prepare_traitrates_run(gene_name, msa_file, char_file, gene_tree_file, species_tree_file, reroot_taxa, base_dir,\
				path_traitrates_exe, path_traitrates_config_template, tree_schema="newick"):
	out_traitrates_dir = base_dir + "/" + gene_name + "/"
	out_msa_file = out_traitrates_dir + gene_name + ".clean.aln"
	out_char_file = out_traitrates_dir + gene_name + ".encoding.txt"
	out_tree_file = out_traitrates_dir + gene_name + ".pruned.tre"
	out_config_file = out_traitrates_dir + gene_name + ".traitrates.cfg"
	
	out_result_file_null = gene_name + ".result.null.txt"
	out_result_file_main = gene_name + ".result.main.txt"
	out_log_file = gene_name + ".log"
	out_ll_pos_file = gene_name + ".ll_pos.txt"
	
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
	
	# output files needed for traitrates run
	system("mkdir " + out_traitrates_dir)
	
	with open(out_tree_file, 'w') as out_fh:
		tree_str = tree_species.as_string(schema=tree_schema).replace("'", '').strip()
		out_fh.write(tree_str)
	
	# prepare encoding file
	with open(char_file, 'r') as in_fh, open(out_char_file, 'w') as out_fh:
		for line in in_fh:
			(taxon, char) = line.strip().split("\t")
			if taxon in taxa_species:
				out_fh.write(">" + taxon + "\n")
				out_fh.write(char + "\n")
	
	# prepare MSA file
	with open(msa_file, 'r') as in_fh, open(out_msa_file, 'w') as out_fh:
		for record in SeqIO.parse(in_fh, 'fasta'):
			clean_record_id = clean_taxon_label(record.id)
			if clean_record_id in taxa_species:
				out_fh.write(">" + clean_record_id + "\n")
				out_fh.write(str(record.seq) + "\n")
	
	# prepare config file
	with open(path_traitrates_config_template, 'r') as in_fh, open(out_config_file, 'w') as out_fh:
		for line in in_fh:
			line = line.strip()
			if line.endswith('TREE_FILE'):
				line = line.replace('TREE_FILE', out_tree_file)
			elif line.endswith('CHAR_FILE'):
				line = line.replace('CHAR_FILE', out_char_file)
			elif line.endswith('MSA_FILE'):
				line = line.replace('MSA_FILE', out_msa_file)
			elif line.endswith('OUT_DIR'):
				line = line.replace('OUT_DIR', out_traitrates_dir)
			elif line.endswith('RESULT_FILE_NULL'):
				line = line.replace('RESULT_FILE_NULL', out_result_file_null)
			elif line.endswith('RESULT_FILE_MAIN'):
				line = line.replace('RESULT_FILE_MAIN', out_result_file_main)
			elif line.endswith('LOG_FILE'):
				line = line.replace('LOG_FILE', out_log_file)
			elif line.endswith('LL_POS_FILE'):
				line = line.replace('LL_POS_FILE', out_ll_pos_file)
			out_fh.write(line + "\n")
	
	traitrates_cmd = " ".join(["cd", out_traitrates_dir, ";", path_traitrates_exe, out_config_file])
	
	return traitrates_cmd


# prepare traitrates files, gene by gene
in_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/"
in_suffix = ".nt.aln.trimmed"
in_species_tree_file = "/data/home/szhan/projects/misc/cp-red/analysis/species_trees_concat/RAxML_bestTree.concat_cds_aa"
in_encoding_file = "/data/home/szhan/projects/misc/cp-red/data/encodings_calcification.txt"
out_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/traitrates/calcification/"
#in_encoding_file = "/data/home/szhan/projects/misc/cp-red/data/encodings_environment.txt"
#out_base_dir = "/data/home/szhan/projects/misc/cp-red/analysis/traitrates/environment/"

traitrates_exe = "/data/home/szhan/projects/misc/cp-red/bin/TraitRateProp_source_20161221/programs/traitRateProp/traitRate.doubleRep"
traitrates_config_template = "/data/home/szhan/projects/misc/cp-red/scripts/traitrates.cfg.template"

in_msa_dir = in_base_dir + "gene_msas/"
in_tree_dir = in_base_dir + "mb_trees/"

reroot_taxa = [	'Cyanidioschyzon_merolae_Cyanidiales_Cyanidiaceae_NC_004799_1',
		'Cyanidium_sp__Cyanidiales_Cyanidiaceae_KJ569775_1',
		'Cyanidium_caldarium_Cyanidiales_Cyanidiaceae_NC_001840_1',
		'Galdieria_partita_Cyanidiales_Galdieriaceae_GPAL',
		'Galdieria_sulphuraria_Cyanidiales_Galdieriaceae_NC_024665_1'
		]

gene_names = [x.replace(in_suffix, '') for x in listdir(in_msa_dir) if x.endswith(in_suffix)]
for in_gene_name in gene_names:
	in_msa_file = in_msa_dir + in_gene_name + in_suffix
	in_gene_tree_file = in_tree_dir + "/" + in_gene_name + "/" + in_gene_name + ".pruned.tre.t"
	traitrates_cmd = prepare_traitrates_run(gene_name=in_gene_name, msa_file=in_msa_file, char_file=in_encoding_file,
						gene_tree_file=in_gene_tree_file, species_tree_file=in_species_tree_file,
						reroot_taxa=reroot_taxa, base_dir=out_base_dir,
						path_traitrates_exe=traitrates_exe, path_traitrates_config_template=traitrates_config_template)
	print traitrates_cmd

