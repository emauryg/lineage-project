## Code to gets shared variants and plot their AAFs.
library(ggplot2)

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
               alt_reads=as.numeric(strsplit(numbers[[1]][2],split=",")[[1]][2])
               ref_reads=as.numeric(numbers[[1]][3])
               AF=alt_reads/ref_reads
               records[[record]] <- c(records[[record]], AF)
          }
     }
     close(vcf_file)
}

count_shared <- lapply(records, function(x) length(x))
shared_vars <- names(count_shared)[which(count_shared>1 & count_shared < length(samples))]
avg_AF <- lapply(records[shared_vars], function(x) mean(x))
AFs<- as.matrix(unlist(avg_AF))

breaks <- pretty(range(AFs), n=nclass.FD(AFs), min.n=1)
bwidth <- breaks[2]-breaks[1]
AFs=as.data.frame(AFs)
pdf("allele_frequency_plot_SNVs_shared.pdf")
ggplot(data=AFs, aes(AFs$V1))+
geom_histogram(col='black', fill='cyan3', alpha=0.5, binwidth=0.01, aes(y=..density..))+ #f for AF40 I used binwidth=0.01
geom_density(col=2)+
xlab("Avg AAF Bulk Phasing Shared")+
ylab("Density")+
theme_bw()
dev.off()
