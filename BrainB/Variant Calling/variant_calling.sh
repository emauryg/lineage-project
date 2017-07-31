## Comand line code of analysis workflow
## To be run blockwise in the command line.
## Make sure to change the extensions depending on the file you want to filter.

module load gcc/6.2.0
module load bcftools/1.3.1
module load annovar/20160201-patched
module load R/3.3.3

#############################################################################################################################
# Find every three phased variant

bulk="/n/data1/bch/genetics/walsh-park/data/Linkage/data_feb_8_2017/1465_1024-pfc-bulk-hc"
chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )

cd ${bulk}
for chr in "${chroms[@]}"
do
     cd ${chr}
     sbatch -p short -n 1 -t 0-1:00:00  --wrap="Rscript /n/scratch2/eam63/merging_cells_project/merging/gvcf_files/three_phased.R ${chr}"
     cd ..
done
## save data in  /n/scratch2/eam63/merging_cells_project/merging/gvcf_files/three_phase_calls_bulkPFC

############################################################################################################################
# Break down the bulk 200x sequence by chromosome

chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )

for chr in "${chroms[@]}"
do
     bcftools view -r ${chr} 1465_1024-pfc-bulk.HC_G.g.vcf.gz -o ${chr}_1465_1024_pfc-bulk.g.vcf
done

############################################################################################################################
# Find SNV candidates from bulk tissue

chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
for i in "${chroms[@]}"
do
     sbatch -p medium -n 1 -t 2-0:00:00 --mem=16G -o ${i}.out  --wrap="Rscript phased_candidates.R ${i}"
done

# Filter out the candidate variants from the gvcf file of bulk. Spread and gather.
chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
for chr in "${chroms[@]}"
do
     bcftools view -R ${chr}_candidate_phased_SNVs.txt 1465_1024-pfc-bulk.HC_G.g.vcf.gz  -o ${chr}_bulk_SNV_candidates.vcf
done

# Gather all SNV candidates
chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
extension="_200X_bulkSNV_candidates.vcf"
cnt=${#chroms[@]}

for ((i=0;i<cnt;i++))
do
     concat_files[i]="${chroms[i]}_bulk_SNV_candidates.vcf"
done

bcftools concat  ${concat_files[@]} -o 1465_ctx${extension}

############################################################################################################################
# Getting the single cell files

files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
extension=".HC_G.g.vcf.gz"
for f in "${files[@]}"
do
     cp ${f}${extension}.tbi /n/scratch2/eam63/merging_cells_project/merging/gvcf_files/${f}/
done

# Breaking down each cell by chromosome
chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
for f in "${files[@]}"
do
     for chr in "${chroms[@]}"
     do
         sbatch -p short -n 1 -t 0-1:00:00 --mem=16G  --wrap="bcftools view -r ${chr} ${f}/${f}.HC_G.g.vcf.gz -o ${f}/${chr}_${f}.g.vcf"
     done
done

############################################################################################################################
# Compare hets from each sinle neurons with the bulk homo refs

chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
for f in "${files[@]}"
do
     for chr in "${chroms[@]}"
     do
          sbatch -p short -n 1 -t 0-6:00:00 --mem=16G  --wrap="Rscript phase_calling.R ${f} ${chr}"
     done
done

# Create filtered out vcf files for each chromosome of each neuron
chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
for f in "${files[@]}"
do
     for chr in "${chroms[@]}"
     do
          sbatch -p short -n 1 -t 0-1:00:00 --mem 16G --wrap="bcftools view -R ${f}/${chr}_${f}_SNVs.txt ${f}/${f}.HC_G.g.vcf.gz  -o ${f}/${chr}_${f}_SNV_candidates.vcf"
     done
done

# Concatenate all the chromosomes so that we have a vcf file per neuron with the SNVs.
chroms=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
files=( "1465-cortex_1-neuron_MDA_3_WGSb" "1465-cortex_1-neuron_MDA_12" "1465-cortex_1-neuron_MDA_2_WGSb" "1465-cortex_1-neuron_MDA_24" "1465-cortex_1-neuron_MDA_6_WGSb" "1465-cortex_1-neuron_MDA_18" "1465-cortex_1-neuron_MDA_39" "1465-cortex_1-neuron_MDA_47" "1465-cortex_1-neuron_MDA_51_WGSb" "1465-cortex_1-neuron_MDA_20"  "1465-cortex_1-neuron_MDA_25" "1465-cortex_1-neuron_MDA_30" "1465-cortex_1-neuron_MDA_43" "1465-cortex_1-neuron_MDA_46" "1465-cortex_1-neuron_MDA_5" "1465-cortex_1-neuron_MDA_8" )
for f in "${files[@]}"
do

     concat_files=${f}/*_${f}_SNV_candidates.vcf

     bcftools concat  ${concat_files[@]} -o ${f}_SingleCell_SNVs.vcf
done


############################################################################################################################



