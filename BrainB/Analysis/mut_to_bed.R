# R code to Convert mutation files into bed files. Notice that the SNV positions are POS-1, POS

samples = c( "1465-cortex_1-neuron_MDA_3_WGSb", "1465-cortex_1-neuron_MDA_12", "1465-cortex_1-neuron_MDA_2_WGSb", "1465-cortex_1-neuron_MDA_24", "1465-cortex_1-neuron_MDA_6_WGSb", "1465-cortex_1-neuron_MDA_18", "1465-cortex_1-neuron_MDA_39", "1465-cortex_1-neuron_MDA_47", "1465-cortex_1-neuron_MDA_51_WGSb", "1465-cortex_1-neuron_MDA_20",  "1465-cortex_1-neuron_MDA_25", "1465-cortex_1-neuron_MDA_30", "1465-cortex_1-neuron_MDA_43", "1465-cortex_1-neuron_MDA_46", "1465-cortex_1-neuron_MDA_5", "1465-cortex_1-neuron_MDA_8" )
extension='_SingleCell_SNVs_1000G_microsat_SegDup_db_noIndels_AF30.vcf'

out_table=c()
for (s in samples){
  vcf_file=file(paste(s,extension,sep=""),open='r')
  lines=readLines(vcf_file)
  for (line in lines){
    if(!grepl('#',line)){
      data=strsplit(line,split="\t")[[1]]
      chrom=data[1]
      pos=as.integer(data[2])
      pos_minus_1=as.character(pos-1)
      ref=data[4]
      alt=gsub(",<NON_REF>","",data[5])
      if (nchar(alt) ==1 & nchar(ref)==1){
        global_id=paste(chrom, pos, ref, alt,sep=";")
        out_table=rbind(out_table,c(paste('chr',chrom,sep=""), pos_minus_1, pos, ref, alt, global_id, s))
      }


    }
  }
  close(vcf_file)
}
write.table(out_table, file=paste("strandbias/ALL_mutations_AF30.bed",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


# R code to create class file
samples = c( "1465-cortex_1-neuron_MDA_3_WGSb", "1465-cortex_1-neuron_MDA_12", "1465-cortex_1-neuron_MDA_2_WGSb", "1465-cortex_1-neuron_MDA_24", "1465-cortex_1-neuron_MDA_6_WGSb", "1465-cortex_1-neuron_MDA_18", "1465-cortex_1-neuron_MDA_39", "1465-cortex_1-neuron_MDA_47", "1465-cortex_1-neuron_MDA_51_WGSb", "1465-cortex_1-neuron_MDA_20",  "1465-cortex_1-neuron_MDA_25", "1465-cortex_1-neuron_MDA_30", "1465-cortex_1-neuron_MDA_43", "1465-cortex_1-neuron_MDA_46", "1465-cortex_1-neuron_MDA_5", "1465-cortex_1-neuron_MDA_8" )

class = "neuron"

record=c()
for (s in samples){
  record = rbind(record, paste(s,class,sep="\t"))
}

write.table(record, file=paste("strandbias/samples_byclass.txt",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



