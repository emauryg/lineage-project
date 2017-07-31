## Comand line code of figures workflow
## To be run blockwise in the command line.
## Make sure to change the extensions depending on the file you want to filter.


module load gcc/6.2.0
module load bcftools/1.3.1
module load annovar/20160201-patched
module load R/3.3.3

#############################################################################################################################
# Plot AAFs of the Single Cell SNVs

Rscript plot_raw_AAFs.R

#############################################################################################################################
# Plot Shared variants AAFs

Rscript plot_shared_AAFs.R

#############################################################################################################################
