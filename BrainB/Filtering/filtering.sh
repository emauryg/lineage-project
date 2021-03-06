## Comand line code of filtering workflow
## To be run blockwise in the command line.
## Make sure to change the extensions depending on the file you want to filter.

module load gcc/6.2.0
module load bcftools/1.3.1
module load annovar/20160201-patched
module load R/3.3.3
module load vcftools/0.1.15
module load htslib/1.3.2

#############################################################################################################################
# Command line code to convert single cell vcf files into annovar input
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
for f in "${files[@]}"
do
     convert2annovar.pl -format vcf4 ${f}_SingleCell_SNVs.vcf -outfile ${f}.avinput
done

#############################################################################################################################
# Identifying the the alleles that have a MAF greater than 0.01 in 1000G database
for f in "${files[@]}"
do
     annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg19 -out ${f} ${f}.avinput humandb -maf 0.01
done

# Get the coordinates to be kept
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension=".hg19_ALL.sites.2015_08_filtered"
out_file_extension="_1000G_maf01.txt"
for f in "${files[@]}"
do
     Rscript get_reads.R ${f} ${extension} ${out_file_extension}
done

# Compress Single Cell vcf files
extension="_SingleCell_SNVs"

for f in "${files[@]}"
do
     vcf-sort -c ${f}${extension}.vcf > ${f}${extension}_sorted.vcf
done

for f in "${files[@]}"
do
      bgzip -c ${f}${extension}_sorted.vcf > ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
      tabix -p vcf -f ${f}${extension}.vcf.gz
done

# Filter out regions of 1000G MAF greater than 0.01
extension="_SingleCell_SNVs"
for f in "${files[@]}"
do
     bcftools view -T ${f}_1000G_maf01.txt ${f}${extension}.vcf.gz -o ${f}${extension}_1000G.vcf
done

#############################################################################################################################
# Filter out microsatellite regions

files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension="_SingleCell_SNVs_1000G.vcf"
out_file_extension="_microsatellites.txt"

for f in "${files[@]}"
do
     sbatch -n 1 -p short -t 0-00:30:00 --mem=16G --wrap="Rscript filter_microsat.R ${f} ${extension} ${out_file_extension}"
done

extension="_SingleCell_SNVs_1000G"

for f in "${files[@]}"
do
     vcf-sort -c ${f}${extension}.vcf > ${f}${extension}_sorted.vcf
done

for f in "${files[@]}"
do
      bgzip -c ${f}${extension}_sorted.vcf > ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
      tabix -p vcf -f ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
     bcftools view -T ${f}_microsatellites.txt ${f}${extension}.vcf.gz -o ${f}${extension}_microsat.vcf
done

#############################################################################################################################
# Filter out segmetnal duplication regions
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension="_SingleCell_SNVs_1000G_microsat.vcf"
out_file_extension="_SegDups.txt"

for f in "${files[@]}"
do
     sbatch -n 1 -p short -t 0-01:00:00 --mem=16G --wrap="Rscript filter_segdup.R ${f} ${extension} ${out_file_extension}"
done

extension="_SingleCell_SNVs_1000G_microsat"

for f in "${files[@]}"
do
     vcf-sort -c ${f}${extension}.vcf > ${f}${extension}_sorted.vcf
done

for f in "${files[@]}"
do
      bgzip -c ${f}${extension}_sorted.vcf > ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
      tabix -p vcf -f ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
     bcftools view -T ${f}${out_file_extension} ${f}${extension}.vcf.gz -o ${f}${extension}_SegDup.vcf
done

#############################################################################################################################
# Filter out duplicated records

moulde load vcflib

files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension="_SingleCell_SNVs_1000G_microsat_SegDup"

for f in "${files[@]}"
do
  vcfuniq ${f}${extension}.vcf > ${f}${extension}_db.vcf
done

#############################################################################################################################
# Remove Indels
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension="_SingleCell_SNVs_1000G_microsat_SegDup_db.vcf"
out_file_extension="_indels.txt"

for f in "${files[@]}"
do
  Rscript filter_indels.R ${f} ${extension} ${out_file_extension}
done

extension="_SingleCell_SNVs_1000G_microsat_SegDup_db"

for f in "${files[@]}"
do
     vcf-sort -c ${f}${extension}.vcf > ${f}${extension}_sorted.vcf
done

for f in "${files[@]}"
do
      bgzip -c ${f}${extension}_sorted.vcf > ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
      tabix -p vcf -f ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
     bcftools view -T ${f}${out_file_extension} ${f}${extension}.vcf.gz -o ${f}${extension}_noIndels.vcf
done

#############################################################################################################################
# Remove AF < 30

files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension="_SingleCell_SNVs_1000G_microsat_SegDup_db_noIndels.vcf"
out_file_extension="_indels.txt"

for f in "${files[@]}"
do
  Rscript filter_AF.R ${f} ${extension} ${out_file_extension}
done

extension="_SingleCell_SNVs_1000G_microsat_SegDup_db_noIndels"

for f in "${files[@]}"
do
     vcf-sort -c ${f}${extension}.vcf > ${f}${extension}_sorted.vcf
done

for f in "${files[@]}"
do
      bgzip -c ${f}${extension}_sorted.vcf > ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
      tabix -p vcf -f ${f}${extension}.vcf.gz
done


for f in "${files[@]}"
do
     bcftools view -T ${f}${out_file_extension} ${f}${extension}.vcf.gz -o ${f}${extension}_AF30.vcf
done






