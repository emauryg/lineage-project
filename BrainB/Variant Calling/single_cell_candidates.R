###############################################################################
### Code to compare hets from each single neuron to the bulk homo refs.
## Output: coordinates of the phased SNVs in each neuron by chromosome

get_bulk_refs <- function( bulk_file_line){
     if (!grepl("#", bulk_file_line)){
          data=strsplit(bulk_file_line, split="\t")[[1]]
          chrom=data[1]
          pos=data[2]
          REF=data[4]
          ALT=gsub(',<NON_REF>', "", data[5])
          record=paste(chrom,pos, REF,ALT, sep=";")
          return(record)
     }
}

find_somatic <- function(bulk_refs, neuron_file_line){
     if (!grepl('#', neuron_file_line) & grepl('0/1', neuron_file_line)){
          data=strsplit(neuron_file_line, split="\t")[[1]]
          chrom=data[1]
          pos=data[2]
          REF=data[4]
          ALT=gsub(',<NON_REF>', "",data[5])
          record=paste(chrom,pos,REF,ALT,sep=";")
          if (record %in% bulk_refs){
               somatic = paste(chrom, pos,sep="\t")
               return(somatic)
          }
     }

}

args <- commandArgs(trailingOnly=TRUE)
neuron_name=args[1]
chr=args[2]
refPath= "/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/three_phase_calls_bulkPFC/"

bulk_file <- file(paste(refPath,chr,'_bulk_SNV_candidates.vcf',sep=""),open='r')
bulk_file_lines <- readLines(bulk_file)

bulk_refs = as.vector(unlist(sapply(bulk_file_lines, function(x) get_bulk_refs(x))))

close(bulk_file)

neuron_file <- file(paste(neuron_name,"/",chr,'_',neuron_name,".g.vcf",sep=""),open='r')
neuron_file_lines <- readLines(neuron_file)

somatics = as.vector(unlist(sapply(neuron_file_lines, USE.NAMES=FALSE, function(x) find_somatic(bulk_refs,x))))
close(neuron_file)

write.table(somatics,file=paste(neuron_name,"/",chr,'_',neuron_name,"_SNVs.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
