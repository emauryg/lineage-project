## Code to make binary SNV matrix with the neuron names on the rows and the shared mutations on the columns

samples = c( "1465-cortex_1-neuron_MDA_3_WGSb", "1465-cortex_1-neuron_MDA_12", "1465-cortex_1-neuron_MDA_2_WGSb", "1465-cortex_1-neuron_MDA_24", "1465-cortex_1-neuron_MDA_6_WGSb", "1465-cortex_1-neuron_MDA_18", "1465-cortex_1-neuron_MDA_39", "1465-cortex_1-neuron_MDA_47", "1465-cortex_1-neuron_MDA_51_WGSb", "1465-cortex_1-neuron_MDA_20",  "1465-cortex_1-neuron_MDA_25", "1465-cortex_1-neuron_MDA_30", "1465-cortex_1-neuron_MDA_43", "1465-cortex_1-neuron_MDA_46", "1465-cortex_1-neuron_MDA_5", "1465-cortex_1-neuron_MDA_8" )
extension='_SingleCell_SNVs.vcf'
records<- vector("list",0)
for (s in samples){
     vcf_file=file(paste(s,extension, sep=""), open='r')
     lines=readLines(vcf_file)
     for (line in lines){
          if(!grepl('#', line)){
               data=strsplit(line, split="\t")[[1]]
               chrom=data[1]
               pos=data[2]
               REF=data[4]
               ALT=gsub(",<NON_REF>","", data[5])
               record=paste(chrom,pos, REF, ALT,sep=";")
               numbers=strsplit(data[10],split=":")
               neuron_name <- gsub("1465-cortex_1-", "", s)
               records[[record]] <- c(records[[record]], neuron_name)
          }
     }
     close(vcf_file)
}

count_shared <- lapply(records, function(x) length(x))
shared_vars <- names(count_shared)[which(count_shared>1 & count_shared < length(samples))]
shared_records <- records[shared_vars]

binary_matrix <- matrix(0, length(samples), length(shared_vars))
rownames(binary_matrix) <- gsub("1465-cortex_1-","", samples)
colnames(binary_matrix) <- shared_vars
for (i in 1:ncol(binary_matrix)){
     mutation <- shared_vars[i]
     cells <- shared_records[[mutation]]
     binary_matrix[cells, i] = 1
}

write.table(binary_matrix, file="binary_matrix_PhaseCalling_BrainB.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
