require(ape)
require(phangorn)


tre_dir <- "../../analysis2/species_trees/"
sp_tree_file_aa <- paste0(tre_dir, "aa_reduced_gamma/RAxML_bestTree.concat_aa_reduced_gamma")
sp_tree_file_nt <- paste0(tre_dir, "nt_reduced/RAxML_bestTree.concat_nt_reduced")

raxml_files <- list.files(path="../../analysis2/gene_trees/aa/", pattern="RAxML_bipartitionsBranchLabels")


tre_dir <- "../../analysis2/gene_trees/"
gene_nrf_mat <- matrix(nrow=length(raxml_files), ncol=1)
gene_names <- c()
for(i in 1:length(raxml_files)){
  gene_name <- gsub("_PROTGAMMA", "", gsub("RAxML_bipartitionsBranchLabels.", "", raxml_files[i]))
  gene_tree_aa <- unroot(read.tree(paste0(tre_dir, "aa/RAxML_bestTree.", gene_name , "_PROTGAMMA")))
  gene_tree_nt <- unroot(read.tree(paste0(tre_dir, "nt/RAxML_bestTree.", gene_name, "_GTRGAMMA")))
  gene_nrf <- RF.dist(gene_tree_aa, gene_tree_nt, normalize=TRUE, check.labels=TRUE, rooted=FALSE)
  gene_names[i] <- gene_name
  gene_nrf_mat[i,1] <- gene_nrf
}
rownames(gene_nrf_mat) <- gene_names


for(i in 1:length(raxml_files)){
 gene_name <- gsub("_PROTGAMMA", "", gsub("RAxML_info.", "", raxml_files[i]))
 gene_tree_aa <- unroot(read.tree(paste0("gene_trees/aa/RAxML_bestTree.", gene_name , "_PROTGAMMA")))
 gene_tree_nt <- unroot(read.tree(paste0("gene_trees/nt/RAxML_bestTree.", gene_name, "_GTRGAMMA")))
 # prune tips in species tree absent in gene tree
 sp_tree_aa <- read.tree(sp_tree_file_aa)
 sp_tree_nt <- read.tree(sp_tree_file_nt)
 tips2prune_aa <- sp_tree_aa$tip.label[!sp_tree_aa$tip.label %in% gene_tree_aa$tip.label]
 tips2prune_nt <- sp_tree_nt$tip.label[!sp_tree_nt$tip.label %in% gene_tree_nt$tip.label]
 sp_tree_aa_pruned <- drop.tip(unroot(sp_tree_aa), tips2prune_aa) # REF
 sp_tree_nt_pruned <- drop.tip(unroot(sp_tree_nt), tips2prune_nt) # REF
 # compute tree distances
 # normalized RF
 nrf_aa_aa <- RF.dist(sp_tree_aa_pruned, gene_tree_aa,
                      normalize=TRUE, check.labels=TRUE, rooted=FALSE)
 nrf_aa_nt <- RF.dist(sp_tree_aa_pruned, gene_tree_nt,
                      normalize=TRUE, check.labels=TRUE, rooted=FALSE)
 nrf_nt_nt <- RF.dist(sp_tree_nt_pruned, gene_tree_nt,
                      normalize=TRUE, check.labels=TRUE, rooted=FALSE)
 nrf_nt_aa <- RF.dist(sp_tree_nt_pruned, gene_tree_aa,
                      normalize=TRUE, check.labels=TRUE, rooted=FALSE)
 # SPR
 spr_aa_aa <- SPR.dist(sp_tree_aa_pruned, gene_tree_aa)
 spr_aa_nt <- SPR.dist(sp_tree_aa_pruned, gene_tree_nt)
 spr_nt_nt <- SPR.dist(sp_tree_nt_pruned, gene_tree_nt)
 spr_nt_aa <- SPR.dist(sp_tree_nt_pruned, gene_tree_aa)
 # print
 print(paste0(c(gene_name,
                nrf_aa_aa, nrf_aa_nt, nrf_nt_nt, nrf_nt_aa,
                spr_aa_aa, spr_aa_nt, spr_nt_nt, spr_nt_aa),
              sep=",", collapse=""))
}


base_path <- "../../analysis2/"
species_tree_path <- paste0(base_path, "species_trees/aa_reduced_gamma/")
gene_tree_path <- paste0(base_path, "gene_trees/aa/")
raxml_prefix <- "RAxML_bootstrap."
raxml_files <- list.files(path=gene_tree_path, pattern=raxml_prefix)

species_tree_file_aa <- paste0(species_tree_path, raxml_prefix, "concat_aa_reduced_gamma")
gene_tree_files_aa <- list.files(path=gene_tree_path, pattern=raxml_prefix)

gene_nrf_matrix <- matrix(nrow=length(raxml_files), ncol=1)
gene_name_list <- c()

for(i in 1:length(raxml_files)){
  gene_name <- gsub("_PROTGAMMA", "", gsub(raxml_prefix, "", raxml_files[i]))
  gene_name_list[i] <- gene_name
  
  gene_tree_aa <- unroot(read.tree(paste0(gene_tree_path, raxml_prefix, gene_name, "_PROTGAMMA")))
  species_tree_aa <- read.tree(species_tree_file_aa)
  
  sample_size <- 100
  random_gene_indices <- sample.int(length(gene_tree_aa), sample_size, replace=F)
  random_species_indices <- sample.int(length(species_tree_aa), sample_size, replace=F)
  
  list_nrf_aa_aa <- c()
  for(j in 1:sample_size){
    random_gene_tree <- gene_tree_aa[[random_gene_indices[j]]]
    random_species_tree <- species_tree_aa[[random_species_indices[j]]]
    # prune tips in species tree absent in gene tree
    tips2prune <- random_species_tree$tip.label[!random_species_tree$tip.label %in% random_gene_tree$tip.label]
    random_species_tree_pruned <- drop.tip(unroot(random_species_tree), tips2prune)
    list_nrf_aa_aa[j] <- RF.dist(random_species_tree_pruned, random_gene_tree,
                                 normalize=TRUE, check.labels=TRUE, rooted=FALSE)
  }
  
  median_nrf_aa_aa <- median(list_nrf_aa_aa)
  gene_nrf_matrix[i,1] <- median_nrf_aa_aa
}

rownames(gene_nrf_matrix) <- gene_name_list

write.table(gene_nrf_matrix, file="test.txt", quote=F, sep="\t")
