require(phytools)

convert_taxon_names <- function(phy) {
 old_info <- phy$tip.label
 new_info <- c()
 for( i in 1:length(old_info) ){
  tax <- unlist(strsplit(old_info[i], split="-"))
  sp <- unlist(strsplit(tax[1], split="_"))
  new_sp <- paste0(sp[1], "_", sp[2], "_", tax[4])
  #new_sp <- paste0("paste(italic('", sp[1], " ", sp[2], "')", ", '", tax[4], "')")
  new_info[i] <- new_sp
 }
 d <- as.data.frame(cbind(label=old_info, nlabel=new_info), stringsAsFactors=F)
 phy$tip.label <- d[[2]][match(phy$tip.label, d[[1]])]
 #phy$tip.label <- sapply(phy$tip.label, function(x) parse(text=x))
 return(phy)
}


tre1 <- read.tree("../analysis2/species_trees/aa_reduced_gamma/RAxML_bestTree.concat_aa_reduced_gamma.tre")
tre2 <- read.tree("../analysis2/gene_trees/aa/RAxML_bestTree.rpoc1_aa2_5end_PROTGAMMA.tre")
#tre2 <- read.tree("../analysis2/gene_trees/aa/RAxML_bestTree.rbcL_PROTGAMMA.tre")

outgroup <- c("Galdieria_sulphuraria-Cyanidiales-Galdieriaceae-NC_024665.1", 
              "Cyanidium_caldarium-Cyanidiales-Cyanidiaceae-NC_001840.1", 
              "Cyanidium_sp.-Cyanidiales-Cyanidiaceae-KJ569775.1", 
              "Cyanidioschyzon_merolae-Cyanidiales-Cyanidiaceae-NC_004799.1")
tre1 <- root(tre1, outgroup)
tre2 <- root(tre2, outgroup)

tre1 <- ladderize(tre1, right=F)
tre2 <- ladderize(tre2, right=F)

phy1 <- convert_taxon_names(tre1)
phy2 <- convert_taxon_names(tre2)
assoc <- cbind(phy1$tip.label, phy1$tip.label)
obj <- cophylo(phy1, phy2, assoc=assoc, print=T, outgroup=T)

pdf("test.pdf")
plot(obj, link.type="curved", link.lwd=1, link.lty="solid", 
     link.col=make.transparent("blue", 0.25), fsize=0.3, pts=0.5)
dev.off()
