library("ggplot2")
library("qqman")

exp <- read.table("mydir/chr/chr21_cnv_counts.txt", header = T, sep = " ", stringsAsFactors = F)

for (i in 20:1)
{
	tmp <- read.table(paste("mydir/chr/chr", i ,"_cnv_counts.txt", sep=""), F, sep = "\t")
	exp <- rbind(tmp, exp)
}

colnames(exp)[1]<-"CHR"
colnames(exp)[2]<-"BP"
colnames(exp)[3]<-"SNP"
colnames(exp)[4]<-"CaseLoss"
colnames(exp)[5]<-"CaseGain"
colnames(exp)[6]<-"ControlLoss"
colnames(exp)[7]<-"ControlGain"
colnames(exp)[8]<-"NoCaseLoss"
colnames(exp)[9]<-"NoControlLoss"
colnames(exp)[10]<-"NoCaseGain"
colnames(exp)[11]<-"NoControlGain"
colnames(exp)[12]<-"percCaseLoss"
colnames(exp)[13]<-"percCaseGain"    
colnames(exp)[14]<-"percControlLoss" 
colnames(exp)[15]<-"percControlGain" 
colnames(exp)[16]<-"moreCaseLosses" 
colnames(exp)[17]<-"rareCaseLoss" 
colnames(exp)[18]<-"moreCaseGains" 
colnames(exp)[19]<-"rareCaseGain" 
colnames(exp)[20]<-"Px.loss"
colnames(exp)[21]<-"Px.gain"

# write to file
# write.table(exp, "mydir/chr/chr_all_cnv_counts.txt" , sep = "\t", row.names = F, col.names = T, quote = F) 

# P loss Chi-square
colnames(exp)[20]<-"P"

png(file="kora-loss.png", width = 900, height = 900)
manhattan(exp, main = "Study Name - Losses", ylim = c(0, 25), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
    chrlabs = c(1:20, "P", "Q"))
dev.off()

png(file="kora-loss-qq.png", width = 350, height = 350)
qq(exp$P, main = "Study Name - Losses Q-Q ", xlim = c(0, 7), ylim = c(0, 
    12), pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()

colnames(exp)[20]<-"Px.loss"

# P gain Chi-square
colnames(exp)[21]<-"P"

png(file="kora-gains.png", width = 900, height = 900)
manhattan(exp, main = "Study Name - Gains", ylim = c(0, 25), cex = 0.6, 
    cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
    chrlabs = c(1:20, "P", "Q"))
dev.off()

png(file="kora-gains-qq.png", width = 350, height = 350)
qq(exp$P, main = "Study Name - Gains Q-Q ", xlim = c(0, 7), ylim = c(0, 
    12), pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()

colnames(exp)[21]<-"Px.gain"
