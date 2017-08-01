# Code to remove the Microsatellite regions


get_microsat_record <- function(lines, mat){
     for (i in 2:length(lines)){
          data=strsplit(lines[i], split="\t")[[1]]
          chrom=gsub("chr","",data[2])
          pos1=as.numeric(data[3])
          pos2=as.numeric(data[4])
          record=c(chrom, pos1, pos2)
          mat  <- rbind(mat, record)


     }
     return(mat)
}

get_single_cell_records <- function(lines, mat){
     for (line in lines){
          if (!grepl('#', line)){
               data=strsplit(line, split="\t")[[1]]
               chrom=data[1]
               pos=as.numeric(data[2])
               record=c(chrom, pos)
               mat <- rbind(mat,record)

          }
     }
     return(mat)

}

get_overlap_position <- function(single_cell, microsat, idx){
     # get the records of the single cell and check if there are interval overlaps
     for (i in 1:nrow(single_cell)){
          pos <- as.numeric(single_cell[i,2])
          chrom <- single_cell[i,1]
          micro <- microsat[which(microsat[,1]==chrom), ]
          for (j in 1:nrow(micro)){
               int1 <- as.numeric(micro[j,2])
               int2 <- as.numeric(micro[j,3])
               if (pos > int1 & pos < int2){
                    idx = rbind(idx, i)
               }
          }
     }
     return(idx)
}



microsat_file <- file("/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/filtering_dbs/microsat_masker.txt.gz", open='r')
lines=readLines(microsat_file)
close(microsat_file)
mat <- matrix(0,0,3)
microsat <- get_microsat_record(lines, mat)

args <- commandArgs(trailingOnly=TRUE)
file_name=args[1]
extension=args[2]
out_file_extension=args[3]

# file_name="1465-cortex_1-neuron_MDA_5"
# extension="_SingleCell_SNVs_1000G_microsat.vcf"


vcf_file <- file(paste(file_name, extension,sep=""), open='r')
lines2=readLines(vcf_file)
close(vcf_file)

mat <- matrix(0,0,2)
single_cell <- get_single_cell_records(lines2, mat)

idx <- c()
filtered_idx <- get_overlap_position(single_cell, microsat, idx)

if (length(filtered_idx) >0){
     filtered_positions <- single_cell[-filtered_idx, ]
     filtered_positions <- unique(filtered_positions)
     write.table(filtered_positions,file=paste(file_name,out_file_extension,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else{
     write.table(single_cell,file=paste(file_name,out_file_extension,sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


