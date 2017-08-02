## Code to filter out indels from the vcf files


get_snvs <- function(line){
     if ( !grepl('#', line) ){
          data=strsplit(line,split="\t")[[1]]
          chrom=data[1]
          pos=data[2]
          REF=data[4]
          ALT=gsub(",<NON_REF>","",data[5])
          if (nchar(ALT) == 1 & nchar(REF)==1){
               record=c(chrom,pos)
               return(record)
          }

     }
}

args <- commandArgs(trailingOnly=TRUE)
file_name=args[1]
extension=args[2]
out_file_extension=args[3]

file_name="1465-cortex_1-neuron_MDA_12"
extension="_SingleCell_SNVs_1000G_microsat_SegDup_db.vcf"



vcf_file <- file(paste(file_name, extension,sep=""),open='r')
lines <- readLines(vcf_file)

filtered_record <- lapply(lines, function(x) get_snvs(x))
filtered_record <- do.call(rbind, filtered_record)

 write.table(filtered_record,file=paste(file_name,out_file_extension,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

