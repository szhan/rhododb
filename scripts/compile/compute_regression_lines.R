result_file <- "../../manuscript/revision2/tableS3.txt"
graph_file <- "Figure_BSnRF-vs-pdist.pdf"


# fit model
dat <- read.table(result_file, head=T, sep="\t")
mod <- lm(bs_nrf_aa_aa_with_replacement ~ median_p_dist_nt, data=dat)

popular_genes <- c("rbcL", "psbA", "psaA", "psaB")
good_genes <- c("gltB", "rpoC1", "rpoB")


pdf(graph_file)
plot(y=dat$bs_nrf_aa_aa_with_replacement,
     x=dat$median_p_dist_nt,
     ylab="Normalized Robinson-Foulds distance (AA trees)",
     xlab="Alignment length (AA)",
     main="", pch=16, col="white")
abline(mod, col="darkgrey", lwd=2)

newx <- seq(0, 1.0, by=0.01)
#newx <- seq(0, 2000, by=100)
pred_interval <- predict(mod, newdata=data.frame(median_p_dist_nt=newx),
                         interval="prediction", level = 0.95)
lines(newx, pred_interval[,2], col="black", lty=2, lwd=2)
lines(newx, pred_interval[,3], col="black", lty=2, lwd=2)

for ( gene in dat$gene ) {
  tmp_pt <- dat[dat$gene == gene,]
  text(y=tmp_pt$bs_nrf_aa_aa_with_replacement,
       x=tmp_pt$median_p_dist_nt, 
       labels=gene, col="grey")
}

for ( gene in popular_genes ) {
  tmp_pt <- dat[dat$gene == gene,]
  text(y=tmp_pt$bs_nrf_aa_aa_with_replacement,
       x=tmp_pt$median_p_dist_nt, 
       labels=gene, col="royalblue")
}

for ( gene in good_genes ) {
  tmp_pt <- dat[dat$gene == gene,]
  text(y=tmp_pt$bs_nrf_aa_aa_with_replacement, 
       x=tmp_pt$median_p_dist_nt, 
       labels=gene, col="coral2")
}

dev.off()
