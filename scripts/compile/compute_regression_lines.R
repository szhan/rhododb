result_file <- "../manuscript/revision2/tableS3.txt"
graph_file <- "Figure4.pdf"


# fit model
dat <- read.table(result_file, head=T, sep="\t")
mod <- lm(nrf_aa_aa ~ median_p_dist_nt, data=dat)

# identify outliers below prediction interval
#pt_pred_int <- predict(mod, newdata=data.frame(median_p_dist_nt=dat$median_p_dist_nt),
#                       interval="prediction", level = 0.95)
#print(dat[dat$nrf_aa_aa < pt_pred_int[,2],]$gene)


popular_genes <- c("rbcL", "psbA", "psaA", "psaB")
good_genes <- c("gltB", "rpoC1", "rpoB")


pdf(graph_file)
plot(y=dat$nrf_aa_aa,
     x=dat$median_p_dist_nt,
     ylab="Normalized Robinson-Foulds distance (AA trees)",
     xlab="p-distance",
     main="", pch=16, col="white")
abline(mod, col="darkgrey", lwd=2)

newx <- seq(0, 1.0, by=0.01)
pred_interval <- predict(mod, newdata=data.frame(median_p_dist_nt=newx),
                         interval="prediction", level = 0.95)
lines(newx, pred_interval[,2], col="black", lty=2, lwd=2)
lines(newx, pred_interval[,3], col="black", lty=2, lwd=2)

for ( gene in dat$gene ) {
  tmp_pt <- dat[dat$gene == gene,]
  text(y=tmp_pt$nrf_aa_aa,
       x=tmp_pt$median_p_dist_nt, 
       labels=gene, col="grey")
}

for ( gene in popular_genes ) {
  tmp_pt <- dat[dat$gene == gene,]
  text(y=tmp_pt$nrf_aa_aa,
       x=tmp_pt$median_p_dist_nt, 
       labels=gene, col="royalblue")
}

for ( gene in good_genes ) {
  tmp_pt <- dat[dat$gene == gene,]
  text(y=tmp_pt$nrf_aa_aa, 
       x=tmp_pt$median_p_dist_nt, 
       labels=gene, col="coral2")
}

dev.off()
