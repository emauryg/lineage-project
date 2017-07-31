Lineage Tracing Project of Normal Brains

The following in the order of the analysis using Phase calling:


# 1)      Identify the 2 SNPs in the site. One will be the germline, “anchor” that was used for 
the phasing analysis.  The other site will be the somatic candidate site.
# 2)      Determine which site is germline and which is a somatic candidate. For a germline 
site to be useful in this analysis, it MUST be heterozygous in the bulk tissue.  We need it to 
have 2 alleles so we can use one of these allele to anchor the somatic SNV on a chromosome. For 
an SNV to be a somatic SNV candidate, it MUST be homozygous in the bulk tissue.  This is 
definitional for us.  We define a somatic SNV as present in a single neurons and ABSENT in the 
bulk. So, the germline anchor must have a “0/1” genotype in the bulk.  You can ignore this 
site.  The somatic candidate site must be “0/” in the bulk.

# 3)      Once you find which site is the somatic SNV site (it is “0/0” in the bulk tissue), 
for every single cell, tally if it is het for this site (“0/1”) or homozygous ref (“0/0”). I 
suppose you would record this information in a matrix with each column corresponding to a 
neuron and each row corresponding to an SNV, and a “1” denoting presence of an SNV, and a “0” 
denoting absence.  We might tweak this in the future.  For example, some single neurons might 
have 15 ref reads and 2 alt reads, and be called “0/0”.  We may want to say that is a hit that 
was just poorly covered.  So keep you code flexible to querying different columns if we want to 
do that later.

# 4)      Transform that matrix into a new one with the neurons both axes, and in each cell 
some measure of co-occurrence of SNVs.

# 5)      Cluster that matrix to find the most related cells
