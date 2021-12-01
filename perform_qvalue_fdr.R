library(qvalue)

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3)  stop("Incorrect number of arguments\n usage\t perform_qvalue_fdr.R INPUT FDR OUTPUT")
opt_input  <- args[1]
opt_fdr    <- as.numeric(args[2])
opt_output <- args[3]

d = read.table(opt_input, hea=FALSE, stringsAsFactors=FALSE)

d=d[!is.na(d$V6),]
d$qval=qvalue(d$V6, lambda=0)$qvalue

write.table(d[which(d$qval <= opt_fdr), ], opt_output, quote=FALSE, row.names=FALSE, col.names=FALSE)
#
