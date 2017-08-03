# Get the number of variants

samples = c( "1465-cortex_1-neuron_MDA_3_WGSb", "1465-cortex_1-neuron_MDA_12", "1465-cortex_1-neuron_MDA_2_WGSb", "1465-cortex_1-neuron_MDA_24", "1465-cortex_1-neuron_MDA_6_WGSb", "1465-cortex_1-neuron_MDA_18", "1465-cortex_1-neuron_MDA_39", "1465-cortex_1-neuron_MDA_47", "1465-cortex_1-neuron_MDA_51_WGSb", "1465-cortex_1-neuron_MDA_20",  "1465-cortex_1-neuron_MDA_25", "1465-cortex_1-neuron_MDA_30", "1465-cortex_1-neuron_MDA_43", "1465-cortex_1-neuron_MDA_46", "1465-cortex_1-neuron_MDA_5", "1465-cortex_1-neuron_MDA_8" )
extension='_SingleCell_SNVs_1000G_microsat_SegDup_db_noIndels_AF30.vcf'

records <-c()
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
               if (!is.element(record,records)){
                    records = rbind(records, record)
               }

          }
     }
     close(vcf_file)
     print(nrow(records))
}
