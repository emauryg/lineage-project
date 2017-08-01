#R code to get the coordinates form a vcf file.

get_records <- function(file_line){
     # this function gets you the coordinates of a mutation in a file
     data=strsplit(file_line, split="\t")[[1]]
     chrom=data[1]
     pos=data[2]
     record=paste(chrom,pos, sep="\t")
     return(record)
}

args <- commandArgs(trailingOnly=TRUE)
file_name=args[1]
extension=args[2]
out_file_extension=args[3]
input_file <- file(paste(file_name,extension,sep=""),open='r')
lines=readLines(input_file)
records <- as.vector(unlist(sapply(lines,USE.NAMES=FALSE, function(x) get_records(x))))
records <- unique(records)
write.table(records,file=paste(file_name,out_file_extension,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
