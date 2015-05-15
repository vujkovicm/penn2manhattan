args <- commandArgs(trailingOnly = TRUE)

# import data
exp.cases <- read.table("/mydir/cases.rawcnv", header = T, sep = " ", stringsAsFactors = F)
exp.controls <- read.table("/mydir/controls.rawcnv", header = T, sep = " ", stringsAsFactors = F)
exp.fam <- read.table("/mydir/plink.fam", header = F, sep = "\t", stringsAsFactors = F)

# remove CNVs with probe length not between 5-100
exp.cases    <- exp.cases[which(exp.cases$size > 4         & exp.cases$size < 101), ]
exp.controls <- exp.controls[which(exp.controls$size > 4   & exp.controls$size < 101), ]

# count # CNV's amd remove peope with >100 CNVS
exp.caseCount    <- as.data.frame(table(exp.cases$regno))
exp.filter.case  <- exp.caseCount[which(exp.caseCount$Freq < 100), ]
exp.exl.case     <- dim(exp.caseCount)[1]    - dim(exp.filter.case)[1]

exp.controlCount <- as.data.frame(table(exp.controls$regno))
exp.filter.control <- exp.controlCount[which(exp.controlCount$Freq < 100), ]
exp.exl.control <- dim(exp.controlCount)[1] - dim(exp.filter.control)[1]

# from fam-file get the total number of patients
exp.nCases    <- table(exp.fam$V6)[[2]]
exp.nControls <- table(exp.fam$V6)[[1]]

# restrict to chr
exp.cases    <- exp.cases[which(exp.cases$regno %in% exp.filter.case$Var1), ]
exp.controls <- exp.controls[which(exp.controls$regno %in% exp.filter.control$Var1), ]

# restrict to chr
exp.chr.cases    <- exp.cases[which(exp.cases$chr == args[[1]]), ]
exp.chr.controls <- exp.controls[which(exp.controls$chr == args[[1]]), ]

# use physical map as output template
out <- read.table(paste("/mydir/chr/chr", args[[1]], ".map", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
colnames(out)[1] <- "CHR" 
colnames(out)[2] <- "BP" 
colnames(out)[3] <- "SNP"

#create output cols
out$exp.GainControl <- out$exp.LossControl <- out$exp.GainCase <- out$exp.LossCase <- 0 

# counting cases
for (i in 1:dim(exp.chr.cases)[1])
{
  if(exp.chr.cases$value[i] < 2 )
  {
    out[which(out$BP == exp.chr.cases$start[i]) : which(out$BP == exp.chr.cases$end[i]), "exp.LossCase"] <- out[which(out$BP == exp.chr.cases$start[i]) : which(out$BP == exp.chr.cases$end[i]), "exp.LossCase"] + 1 
  }
  if(exp.chr.cases$value[i] > 2 )
  {
    out[which(out$BP == exp.chr.cases$start[i]) : which(out$BP == exp.chr.cases$end[i]), "exp.GainCase"] <- out[which(out$BP == exp.chr.cases$start[i]) : which(out$BP == exp.chr.cases$end[i]), "exp.GainCase"] + 1 
  }
}

# counting controls
for (i in 1:dim(exp.chr.controls)[1])
{
  if(exp.chr.controls$value[i] < 2 )
  {
    out[which(out$BP == exp.chr.controls$start[i]) : which(out$BP == exp.chr.controls$end[i]), "exp.LossControl"] <- out[which(out$BP == exp.chr.controls$start[i]) : which(out$BP == exp.chr.controls$end[i]), "exp.LossControl"] + 1 
  }
    if(exp.chr.controls$value[i] > 2 )
  {
    out[which(out$BP == exp.chr.controls$start[i]) : which(out$BP == exp.chr.controls$end[i]), "exp.GainControl"] <- out[which(out$BP == exp.chr.controls$start[i]) : which(out$BP == exp.chr.controls$end[i]), "exp.GainControl"] + 1 
  }
}

# counting counterfactuals and summaries
for (i in 1:dim(out)[1])
{
  out$exp.NoLossCase[i]      <- exp.nCases    - out$exp.LossCase[i]    - exp.exl.case    
  out$exp.NoLossControl[i]   <- exp.nControls - out$exp.LossControl[i] - exp.exl.control 
  out$exp.NoGainCase[i]      <- exp.nCases    - out$exp.GainCase[i]    - exp.exl.case 
  out$exp.NoGainControl[i]   <- exp.nControls - out$exp.GainControl[i] - exp.exl.control  
}

# create percentages
out$exp.percCaseLoss    <- out$exp.LossCase    / (out$exp.LossCase    + out$exp.NoLossCase)
out$exp.percCaseGain    <- out$exp.GainCase    / (out$exp.GainCase    + out$exp.NoGainCase)
out$exp.percControlLoss <- out$exp.LossControl / (out$exp.LossControl + out$exp.NoLossControl)
out$exp.percControlGain <- out$exp.GainControl / (out$exp.GainControl + out$exp.NoGainControl)

# identify rare losses for downstream analysis (<3%)
out$exp.rareCaseGain <- out$exp.moreCaseGains <- out$exp.rareCaseLoss <- out$exp.moreCaseLosses <- 0 

for (i in 1:dim(out)[1])
{
  # losses cases > controls
  if (out$exp.percCaseLoss[i] > out$exp.percControlLoss[i])
  {
    out$exp.moreCaseLosses[i] <- 1
  }
    # losses cases < 1.5%
  if (out$exp.percCaseLoss[i] < 0.03)
  {
    if(out$exp.percCaseLoss[i] > 0)
    {
      out$exp.rareCaseLoss[i] <- 1
    }
  }
  # gains cases > controls
  if (out$exp.percCaseGain[i] > out$exp.percControlGain[i])
  {
    out$exp.moreCaseGains[i] <- 1
  }
  # gains cases < 1.5%
  if (out$exp.percCaseGain[i] < 0.03)
  {
    if(out$exp.percCaseGain[i] > 0 )
    {
      out$exp.rareCaseGain[i] <- 1
    }
  }
}

# Chi-X + Fisher's Exact
out$exp.P.loss <- 9
out$exp.P.gain <- 9

for (i in 1:dim(out)[1])
{
  # all p-values
  out$exp.P.loss[i] <-  chisq.test(rbind(c(out$exp.LossCase[i], out$exp.NoLossCase[i]), c(out$exp.LossControl[i], out$exp.NoLossControl[i])))$p.value 
  out$val.P.loss[i] <-  chisq.test(rbind(c(out$val.LossCase[i], out$val.NoLossCase[i]), c(out$val.LossControl[i], out$val.NoLossControl[i])))$p.value 
}

out$exp.P.loss[out$exp.P.loss == 9] <- NA
out$exp.P.gain[out$exp.P.gain == 9] <- NA

# Manhattan plot
# png(file = paste("/mydir/out/chr", args[[1]], "_cnv_manhattan.png", sep = ""), width = 900, height = 350)
# manhattan(out, colors = c("black", "#666666", "#CC6600"), pch = 20, main = paste ("Chr ", args[[1]], sep = ""), genomewideline = F, suggestiveline = F)
# dev.off()

# QQ-plot
# png(file = paste("/mydir/out/chr", args[[1]], "_cnv_qq.png", sep = ""), width = 300, height = 300)
# qq(out$P)
# dev.off()

# Export
write.table(out, paste("mydir/chr/chr", args[[1]], "_cnv_counts.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F) 

