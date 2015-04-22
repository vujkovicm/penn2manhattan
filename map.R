# Transform map file into chromosome specific files

# import map (plink)
map <- read.table("/mydir/mystudy.map", sep = "\t", F, stringsAsFactors = F)

# keep chromosome and basepair position
map <- map[ , c(1, 4)]

# split into chromosomes and export
for (i in 1:22)
{
  chr <- map[which(map[1] == i, ]
  chr <- chr[order(chr[2]), ] 
  write.table(chr, paste("mydir/chr/chr", i, ".map" sep = ""), row.names = F, col.names = F, quote = F, sep = " ")
}
