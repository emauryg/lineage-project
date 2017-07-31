##########################################################
# Finding every three phase call
########################################################

### working one chromosome at a time

load('linkage.rda')
args<- commandArgs(trailingOnly=TRUE)
chr <- args[1]
path <- '/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/three_phase_calls_bulkPFC/'
load("linkage.rda")
linkage_copy <- linkage
filtered_linkage<- matrix(0,0,7)
colnames(filtered_linkage) <- colnames(linkage)

for (i in 1:nrow(linkage_copy)){
     if (linkage_copy[i,1] >0 & linkage_copy[i,2]>0 & linkage_copy[i,3]>0){
          filtered_linkage = rbind(filtered_linkage, linkage_copy[i,])
          print(paste('row ', i, ' of ', nrow(linkage_copy)))
     }
}

save(filtered_linkage, file=paste(path,"filtered_linkage_chr",chr,".rda",sep=""))
