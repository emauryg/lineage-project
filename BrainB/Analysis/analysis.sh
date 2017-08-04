## Comand line code of analysis workflow
## To be run blockwise in the command line.
## Make sure to change the extensions depending on the file you want to filter.

module load gcc/6.2.0
module load bcftools/1.3.1
module load annovar/20160201-patched
module load R/3.3.3
module load vcftools


#############################################################################################################################
# Make binarry matrix with the neuron names on the rows and the shared mutations on the columns

Rscript make_binary_matrix.R

#############################################################################################################################
# Make binarry matrix with the neuron names on the rows and the shared mutations on the columns after 1000G filtering

Rscript make_binary_matrix_1000G.R

#############################################################################################################################
# SNS count

Rscript mut_to_bed.R

python sns_count.py

#############################################################################################################################
# Strandbias

 bedtools intersect -a ALL_mutations_AF30.bed -b /n/scratch2/eam63/Schizophrenia/mutect/annovar_processing/UCSC.hg19.RefGene.downloaded20171907.bed -wa -wb > ALL_mutations_formatted_AF30.bed

python strandbias.pyt

#############################################################################################################################
# Annovar annotation
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )


extension="_SingleCell_SNVs_1000G_microsat_SegDup_db_noIndels_AF30"
for f in "${files[@]}"
do
     convert2annovar.pl -format vcf4 ${f}${extension}.vcf -outfile ${f}${extension}.avinput
done

for f in "${files[@]}"
do
     annotate_variation.pl --geneanno -dbtype refGene -out  ${f}${extension} -build hg19 ${f}${extension}.avinput  humandb/
done
