# Code to remove the Microsatellite regions

# chrom =data[2]
# pos1=data[3]
# pos2= data[4]

get_microsat_record <- function(line){
     data=strsplit(line, split="\t")[[1]]
     chrom=gsub("chr","",data[2])
     pos1=data[3]
     pos2=data[4]
     record = c(chrom, pos1, pos2)
}

create_microsat_table <- function(lines){
     records <- sapply(lines[2:length(lines)], USE.NAMES=FALSE, function(x) get_microsat_record(x))
}

microsat_file <- file("/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/filtering_dbs/microsat_masker.txt.gz", open='r')
lines=readLines(microsat_file)
micro_records <- create_microsat_table(lines)
