args <- commandArgs(trailingOnly = TRUE)

# import data
cases <- read.table("/mydir/cases.rawcnv", header = T, sep = " ", stringsAsFactors = F)
controls <- read.table("/mydir/controls.rawcnv", header = T, sep = " ", stringsAsFactors = F)
fam <- read.table("/mydir/plink.fam", header = F, sep = "\t", stringsAsFactors = F)
out <- read.table(paste("/mydir/chr/chr", args[[1]], ".map", sep = ""), header = T, sep = "\t", stringsAsFactors = F)

# restrict to chr
cases <- cases[which(cases$chr == args[[1]]), ]
controls <- controls[which(controls$chr == args[[1]]), ]

# use physical map as output template
colnames(out)[1] <- "CHR" 
colnames(out)[2] <- "BP" 
colnames(out)[3] <- "SNP"

# create output cols
out$Xp.loss <- 0 ;
out$Xp.gain <- 0 ;
out$Fp.loss <- 0 ;
out$Fp.gain <- 0 ;
out$LossCase <- 0 ;
out$LossControl <- 0 ;
out$GainCase <- 0 ;
out$GainControl <- 0 ;

# counting cases
for (i in 1:dim(cases)[1])
{
  # first select the correct chromosome for internal count revision
  if(cases$value[i] < 2 )
  {
    out[which(out$BP == cases$start[i]) : which(out$BP == cases$end[i]), "LossCase"] <- out[which(out$BP == cases$start[i]) : which(out$BP == cases$end[i]), "LossCase"] + 1 ;
  }
  if(cases$value[i] > 2 )
  {
    out[which(out$BP == cases$start[i]) : which(out$BP == cases$end[i]), "GainCase"] <- out[which(out$BP == cases$start[i]) : which(out$BP == cases$end[i]), "GainCase"] + 1 ;
  }
  # save new internal count to external data frame
}

# counting controls
for (i in 1:dim(controls)[1])
{
  if(controls$value[i] < 2 )
  {
    out[which(out$BP == controls$start[i]) : which(out$BP == controls$end[i]), "LossControl"] <- out[which(out$BP == controls$start[i]) : which(out$BP == controls$end[i]), "LossControl"] + 1 ;
  }
  if(controls$value[i] > 2 )
  {
    out[which(out$BP == controls$start[i]) : which(out$BP == controls$end[i]), "GainControl"] <- out[which(out$BP == controls$start[i]) : which(out$BP == controls$end[i]), "GainControl"] + 1 ;
  }
}

# Chi-X + Fisher's Exact
for (j in 1:dim(out)[1])
{
  out$Xp.loss[j] <-  chisq.test(rbind(c(out$LossCase[j], table(fam$V6)[[2]] - out$LossCase[j]), c(out$LossControl[j], table(fam$V6)[[1]] - out$LossControl[j])))$p.value ;
  out$Xp.gain[j] <-  chisq.test(rbind(c(out$GainCase[j], table(fam$V6)[[2]] - out$GainCase[j]), c(out$GainControl[j], table(fam$V6)[[1]] - out$GainControl[j])))$p.value ;
  out$Fp.loss[j] <- fisher.test(rbind(c(out$LossCase[j], table(fam$V6)[[2]] - out$LossCase[j]), c(out$LossControl[j], table(fam$V6)[[1]] - out$LossControl[j])))$p.value ;
  out$Fp.gain[j] <- fisher.test(rbind(c(out$GainCase[j], table(fam$V6)[[2]] - out$GainCase[j]), c(out$GainControl[j], table(fam$V6)[[1]] - out$GainControl[j])))$p.value ;
}

# Manhattan plot
# png(file = paste("/mydir/out/chr", args[[1]], "_cnv_manhattan.png", sep = ""), width = 900, height = 350)
# manhattan(out, colors = c("black", "#666666", "#CC6600"), pch = 20, main = paste ("Chr ", args[[1]], sep = ""), genomewideline = F, suggestiveline = F)
# dev.off()

# QQ-plot
# png(file = paste("/mydir/out/chr", args[[1]], "_cnv_qq.png", sep = ""), width = 300, height = 300)
# qq(out$P)
# dev.off()

# Export
write.table(out, paste("mydir/chr/chr", args[[1]], "_cnv_counts.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F) ;

