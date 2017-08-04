## Variant calling command line scripts
## Please run block wise

module load gcc/6.2.0
module load bcftools/1.3.1
module load annovar/20160201-patched
module load R/3.3.3
module load vcftools
module load gatk/3.7

sbatch  -n 4 -p medium -t 2-00:00:00  --mem=16G --wrap="java -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /groups/walsh/genomes/human_g1k_v37_decoy.fasta --variant 1465-cortex_1-neuron_MDA_12.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_18.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_20.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_24.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_25.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_2_WGSb.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_30.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_39.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_3_WGSb.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_43.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_46.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_47.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_51_WGSb.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_5.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_6_WGSb.HC_G.g.vcf.gz --variant 1465-cortex_1-neuron_MDA_8.HC_G.g.vcf.gz --variant 1465_1024-pfc-bulk.HC_G.g.vcf.gz -o all_neurons_variants.g.vcf.gz -nct 4"
