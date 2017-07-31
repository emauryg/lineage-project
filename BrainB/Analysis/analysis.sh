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
