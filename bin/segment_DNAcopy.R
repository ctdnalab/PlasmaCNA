library(DNAcopy)
library(optparse)

option_list <- list(
  make_option(c("--IN"), type = "character", help = "Path to input table. Required."),
  make_option(c("--OUT"), type = "character", help = "Path to output file. Required.")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
options(scipen=0, stringsAsFactors=FALSE)

inTable <- opt$IN
outFile <- opt$OUT

df <-read.table(inTable,sep='\t',header=TRUE)
cna <- CNA(df$depth,df$chrom,df$startPos,presorted=TRUE,data.type="logratio")
cna <- smooth.CNA(cna)
s <- segment(cna)
write.table(s[2],outFile,sep='\t',quote=FALSE,row.names=FALSE)