require(ape)
require(phangorn)


sp_tree_file_aa <- "../analysis2/species_trees/aa_reduced_gamma/RAxML_bestTree.concat_aa_reduced_gamma.tre"
sp_tree_file_nt <- "../analysis2/species_trees/nt_reduced/RAxML_bestTree.concat_nt_reduced"

raxml_files <- list.files(path="gene_trees/aa/", pattern="RAxML_info.")
#raxml_files <- raxml_files[1]

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


# subset tree analysis
raxml_files <- list.files(path="../analysis2/subsets/aa/", pattern="RAxML_bestTree.")
#raxml_files <- raxml_files[1]

for(i in 1:length(raxml_files)){
  gene_name <- gsub(".subset.aa", "", gsub("RAxML_bestTree.", "", raxml_files[i]))
  gene_tree_aa <- unroot(read.tree(paste0("../analysis2/subsets/aa/RAxML_bestTree.", gene_name, ".subset.aa")))
  gene_tree_nt <- unroot(read.tree(paste0("../analysis2/subsets/nt/RAxML_bestTree.", gene_name, ".subset.nt")))
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
  nrf_nt_nt <- RF.dist(sp_tree_nt_pruned, gene_tree_nt,
                       normalize=TRUE, check.labels=TRUE, rooted=FALSE)
  print(paste0(c(gene_name,
                 nrf_aa_aa, nrf_nt_nt),
               sep=",", collapse=""))
}


reduced_set <- c('apcE', 'gltB', 'rpoB', 'rpoc1', 'rpoc2')
#pos_tre <- unroot(read.tree("../analysis2/gene_trees/nt/RAxML_bestTree.rpoc2_GTRGAMMA"))
pos_tre <- unroot(read.tree("../analysis2/min_set/aa/RAxML_bestTree.min_set_aa"))
sp_tre <- unroot(read.tree(sp_tree_file_aa))
tips2prune <- sp_tre$tip.label[!sp_tre$tip.label %in% pos_tre$tip.label]
sp_tre_pruned <- drop.tip(sp_tre, tips2prune)
nrf <- RF.dist(sp_tre_pruned, pos_tre, normalize=TRUE, check.labels=TRUE, rooted=FALSE)
nrf


# get sliding window results
rpoC1_nt <- seq(0, 1170, 90)
rpoC1_aa <- seq(0, 390, 30)
rpoB_nt <- seq(0, 2700, 90)
rpoB_aa <- seq(0, 900, 30)


get_nrf_profile <- function(pos_list, in_prefix, sp_tre_file){
  nrf_list <- c()
  for(i in 1:length(pos_list)){
    pos_tre_file <- paste0(in_prefix, pos_list[i])
    pos_tre <- unroot(read.tree(pos_tre_file))
    sp_tre <- unroot(read.tree(sp_tre_file))
    tips2prune <- sp_tre$tip.label[!sp_tre$tip.label %in% pos_tre$tip.label]
    sp_tre_pruned <- drop.tip(sp_tre, tips2prune)
    nrf <- RF.dist(sp_tre_pruned, pos_tre, normalize=TRUE, check.labels=TRUE, rooted=FALSE)
    nrf_list[i] <- nrf
  }
  return(nrf_list)
}

rpoC1_nt_nrf_window <- get_nrf_profile(rpoC1_nt, "../analysis2/primer_design/RAxML_bestTree.rpoC1.nt.", sp_tree_file_nt)
rpoC1_aa_nrf_window <- get_nrf_profile(rpoC1_aa, "../analysis2/primer_design/RAxML_bestTree.rpoC1.aa.", sp_tree_file_aa)
rpoC1_nt_nrf_ref <- 0.2079208
rpoC1_aa_nrf_ref <- 0.2376238

pdf("FigureS3_rpoC1.pdf")
plot(x=rpoC1_nt, y=rpoC1_nt_nrf_window,
     xlab="Position (NT)", ylab="Normalized Robinson-Foulds distance",
     col="coral2", type="l", ylim=c(0,1))
points(x=rpoC1_nt, y=rpoC1_aa_nrf_window, col="royalblue", type="l")
abline(h=rpoC1_nt_nrf_ref, col="coral2", lty=2)
abline(h=rpoC1_aa_nrf_ref, col="royalblue", lty=2)
dev.off()


rpoB_nt_nrf_window <- get_nrf_profile(rpoB_nt, "../analysis2/primer_design/RAxML_bestTree.rpoB.nt.", sp_tree_file_nt)
rpoB_aa_nrf_window <- get_nrf_profile(rpoB_aa, "../analysis2/primer_design/RAxML_bestTree.rpoB.aa.", sp_tree_file_aa)
rpoB_nt_nrf_ref <- 0.1553398
rpoB_aa_nrf_ref <- 0.2038835

pdf("FigureS4_rpoB.pdf")
plot(x=rpoB_nt, y=rpoB_nt_nrf_window,
     xlab="Position (NT)", ylab="Normalized Robinson-Foulds distance",
     col="coral2", type="l", ylim=c(0,1))
points(x=rpoB_nt, y=rpoB_aa_nrf_window, col="royalblue", type="l")
abline(h=rpoB_nt_nrf_ref, col="coral2", lty=2)
abline(h=rpoB_aa_nrf_ref, col="royalblue", lty=2)
dev.off()




pos_tre <- unroot(read.tree("../analysis2/gene_trees/nt/RAxML_bestTree.rpoB_GTRGAMMA"))
sp_tre <- unroot(read.tree(sp_tree_file_nt))
tips2prune <- sp_tre$tip.label[!sp_tre$tip.label %in% pos_tre$tip.label]
sp_tre_pruned <- drop.tip(sp_tre, tips2prune)
nrf <- RF.dist(sp_tre_pruned, pos_tre, normalize=TRUE, check.labels=TRUE, rooted=FALSE)
