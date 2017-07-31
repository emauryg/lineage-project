###############################################
### R script to find SNV candidates from bulk tissue
###############################################

## SNV candidates are homozygous reference in the bulk tissue. The anchor mutation is heterozygous in the bulk as it is probably germline.


get_homo_refs <- function(line){
     ## This function checks a line of a vcf file
     ## outpus the chrom, pos homo_refs

     if (grepl("0/0", line) & !grepl("#", line)){
          data=strsplit(line, split="\t")[[1]]
          chrom=data[1]
          pos=data[2]
          record=paste(chrom, pos)
          return(record)
     }
}

get_hets <- function(line){
     ## This function checks a line of a vcf file
     ## Outpus the chrom, pos of hets

     if (grepl("0/1",line) & !grepl("#", line)){
          data=strsplit(line, split="\t")[[1]]
          chrom=data[1]
          pos=data[2]
          record=paste(chrom,pos)
          return (record)
     }
}


get_SNV_candidates <- function(hets, homo_refs, filtered_linkage){

     ## This function takes in the output of the get_homo_refs and get_hets functions along with the linkage information
     ##  Outputs the positions of the SNVs candidates and the phased pairs as text files.

     SNVs=c()
     SNPs=c()
     SNV_records=c()

     for (i in 1:nrow(filtered_linkage)){
          data=filtered_linkage[i,]

          site1=as.character(data[5])
          chrom1=strsplit(site1,split=";")[[1]][1]
          pos1=strsplit(site1,split=";")[[1]][2]
          record1=paste(chrom1,pos1)

          site2=as.character(data[6])
          chrom2=strsplit(site2,split=";")[[1]][1]
          pos2=strsplit(site2, split=";")[[1]][2]
          record2=paste(chrom2,pos2)

          if ( is.element(record1, homo_refs) & is.element(record2,hets) ){
               SNVs=rbind(SNVs,site1)
               SNV_records=rbind(SNV_records,c(chrom1, pos1))
               SNPs=rbind(SNPs,site2)
               print(site1)
          } else if( is.element(record2, homo_refs) & is.element(record1,hets)){
               SNVs=rbind(SNVs,site2)
               SNV_records=rbind(SNV_records,c(chrom2,pos2))
               SNPs=rbind(SNPs,site1)
               print(site2)
          }


     }
     phased=cbind(SNVs,SNPs)
     colnames(phased)=c("SNVs", "SNPs")

     write.table(SNV_records,file=paste(chr, "_candidate_phased_SNVs.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
     write.table(phased, file=paste(chr,"phased_pairs.txt",sep=""),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

}

## Runs the above function in order
args<- commandArgs(trailingOnly=TRUE)
chr=args[1]
load(paste("filtered_linkage_chr",chr,".rda",sep=""))
vcf_file_path= '_1465_1024_pfc-bulk.g.vcf'

vcf_file=file(paste(chr, vcf_file_path,sep=""), open='r')
lines=readLines(vcf_file)
homo_refs=as.vector(unlist(sapply(lines, function(x) get_homo_refs(x))))
hets=as.vector(unlist(sapply(lines, function(x) get_hets(x))))

close(vcf_file)
get_SNV_candidates(hets, homo_refs, filtered_linkage)
