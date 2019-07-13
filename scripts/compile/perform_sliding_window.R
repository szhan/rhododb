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
