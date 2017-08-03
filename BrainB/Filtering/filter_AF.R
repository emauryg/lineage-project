# Remove vars that have an allele fraction less 30.

get_AFs <- function(lines){
     records= c()
     for (line in lines){
          if(!grepl('#', line)){
               data=strsplit(line, split="\t")[[1]]
               chrom=data[1]
               pos=data[2]
               REF=data[4]
               ALT=gsub(",<NON_REF>","", data[5])
               record=paste(chrom,pos, sep="\t")
               numbers=strsplit(data[10],split=":")
               alt_reads=as.numeric(strsplit(numbers[[1]][2],split=",")[[1]][2])
               ref_reads=as.numeric(numbers[[1]][3])
               AF=alt_reads/ref_reads
               if (AF > 0.30){
                    records = rbind(records, record)
               }

          }
     }
     return(records)

}

args <- commandArgs(trailingOnly=TRUE)
file_name=args[1]
extension=args[2]
out_file_extension=args[3]

# file_name="1465-cortex_1-neuron_MDA_12"
# extension="_SingleCell_SNVs_1000G_microsat_SegDup_db_noIndels.vcf"

vcf_file=file(paste(file_name,extension, sep=""), open='r')
lines=readLines(vcf_file)

filtered_record <- get_AFs(lines)

 write.table(filtered_record,file=paste(file_name,out_file_extension,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

