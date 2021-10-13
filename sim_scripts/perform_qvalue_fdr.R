library(qvalue)

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4)  stop("Incorrect number of arguments\n usage\t perform_qvalue_fdr.R INPUT FDR METHOD OUTPUT")
opt_input  <- args[1]
opt_fdr    <- as.numeric(args[2])
opt_method <- args[3]
opt_output <- args[4]

#d = read.table("results.genes.full.txt.gz", hea=FALSE, stringsAsFactors=FALSE)
#d=d[!is.na(d$V19),]
#d$qval=qvalue(d$V19, lambda=0)$qvalue
#write.table(d[which(d$qval <= 0.05), ], "results.genes.significant.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

d = read.table(opt_input, hea=FALSE, stringsAsFactors=FALSE)

if(opt_method == 'IsoQTL' || opt_method == 'fisher_perm') {
        d=d[!is.na(d$V6),]
        d$qval=qvalue(d$V6, lambda=0)$qvalue
} else if(opt_method == 'QTLtools') {
        d=d[!is.na(d$V20),]
        d$qval=qvalue(d$V20, lambda=0)$qvalue
} else if(opt_method == 'combined_qtltools') {
        d=d[!is.na(d$V2),]
        d$qval=qvalue(d$V2, lambda=0)$qvalue
}

write.table(d[which(d$qval <= opt_fdr), ], opt_output, quote=FALSE, row.names=FALSE, col.names=FALSE)
#
