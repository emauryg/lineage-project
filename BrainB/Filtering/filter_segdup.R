## R-code to filter out Segemental duplication regions


setwd("/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/")

get_dups_record <- function(line, vcf_records){
     if( !grepl('#', line) ){
          data=strsplit(line, split="\t")[[1]]
          chrom=gsub('chr','',data[2])
          pos1=as.numeric(data[3])
          pos2=as.numeric(data[4])
          #record=c(chrom, pos1, pos2)
          records = invisible(apply(vcf_records,1, function(x) check_interval(chrom, pos1, pos2, x)))

          return(records)
     }
}

check_interval <- function(chrom, pos1, pos2, row){
     chr = row[1]
     pos = as.numeric(row[2])
     pos1=as.numeric(pos1)
     pos2=as.numeric(pos2)

     if ((chr == chrom) & (pos > pos1 & pos < pos2)){
          record = c(chr,pos)
     }

}
get_vcf_record <- function(line){
     if (!grepl('#', line)){
          data=strsplit(line, split="\t")[[1]]
          chrom=data[1]
          pos=data[2]
          record= c(chrom,pos)
     }
}

args <- commandArgs(trailingOnly=TRUE)
file_name=args[1]
extension=args[2]
out_file_extension=args[3]

# file_name="1465-cortex_1-neuron_MDA_43"
# extension="_SingleCell_SNVs_1000G_microsat.vcf"

vcf_file <- file(paste(file_name, extension,sep=""), open='r')
lines <- readLines(vcf_file)

vcf_records <- lapply(lines, function(x) get_vcf_record(x))
vcf_records <- do.call(rbind, vcf_records)

dup_path <- "/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/filtering_dbs/segmental_dups.txt"
dup_file <- file(dup_path, open='r')
lines <- readLines(dup_file)

dups_records <- lapply(lines, function(x) get_dups_record(x, vcf_records))
dups_records <- do.call(rbind, do.call(rbind, dups_records))
dups_records <- unique(dups_records)

records=c()
for (i in 1:nrow(dups_records)){
     chrom=dups_records[i,1]
     pos=dups_records[i,2]
     record=paste(chrom,pos)
     records=rbind(records, record)
}

filtered_segdups =c()
for (i in 1:nrow(vcf_records)){
     chrom=vcf_records[i,1]
     pos=vcf_records[i,2]
     record=paste(chrom, pos)
     if (!is.element(record, records)){
          filtered_segdups=rbind(filtered_segdups, paste(chrom, pos,sep="\t"))
     }
}

 write.table(filtered_segdups,file=paste(file_name,out_file_extension,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


