######################################################################################################################
##
## Python code to analyze Single Nucleotide Substitutions
##
######################################################################################################################

import os
import glob
import numpy as np

def count_sub_types(fname, t_outfname, f_outfname):

    ''' (file)

    Counts the different types of SNVs in a VCF file, so the number of
    C>T, C>A, G>T, etc. can be compared. Makes an outfile with total counts
    and fraction of counts within a cell.
    '''

    file_in = open (fname)
    t_outfile = open(t_outfname, 'w')
    f_outfile = open(f_outfname, 'w')
    tab = '\t'

    # Create dictionaries to count different SNV subtypes

    sub_totals = {'C>A':0, 'C>G':0, 'C>T':0, 'T>A':0, 'T>C':0, 'T>G':0}

    sub_totals_list = sorted(sub_totals)

    rev_comps = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    cell_classes = {}

    sub_totals_by_cell = {}
    totals_by_cell = {}

    # Make a header
    for r in sub_totals_list:
        t_outfile.write (tab + r)
        f_outfile.write (tab + r)
    t_outfile.write ('\n')
    f_outfile.write ('\n')

    for line in file_in:

        # Tally the mutations in each cell in the totals_by_cell
        # and sub_totals_by_cell dictionaries

        if 'chromosome' not in line:
            sline = line.split('\t')
            cell_ID = sline[-1].strip('\n')

            if cell_ID in sub_totals_by_cell:
                totals_by_cell[cell_ID] +=1
            else:
                sub_totals_by_cell[cell_ID] = {'C>A':0, 'C>G':0, 'C>T':0, 'T>A':0, 'T>C':0, 'T>G':0}
                totals_by_cell[cell_ID] = 1


            REF = sline[3]
            ALT = sline[4]

            if REF != 'C' and REF != 'T':
                REF_comp = rev_comps[REF]
                ALT_comp = rev_comps[ALT]
                sub = REF_comp + '>' + ALT_comp

            else:
                sub = REF + '>' + ALT

            sub_totals_by_cell[cell_ID][sub] += 1



    sorted_cell_IDs = sorted(sub_totals_by_cell)

    # Write the outfiles using the data from totals_by_cell and
    # sub_totals_by_cell

    for c in sorted_cell_IDs:
        t_outfile.write(c)
        f_outfile.write(c)
        for s in sub_totals_list:
            t_outfile.write (tab + str(sub_totals_by_cell[c][s]))
            f_outfile.write (tab + str(float(sub_totals_by_cell[c][s])/float(totals_by_cell[c])))

        t_outfile.write('\n')
        f_outfile.write('\n')

    file_in.close()
    t_outfile.close()
    f_outfile.close()

def SNS_grouper(counted_fname, fraction_fname, class_key_fname, info_column,
                tSNS_by_class_outfname, fmsdSNS_by_class_outfname):
    """
    Group quantified data by class, such as age group
    tissue, and disease. Calculate mean and SD for fraction
    data per group.
    """

    fraction_file = open(fraction_fname)
    counted_file = open(counted_fname)
    class_key_file = open(class_key_fname)
    tSNS_by_class_outfile = open(tSNS_by_class_outfname, 'w')
    fmsdSNS_by_class_outfile = open(fmsdSNS_by_class_outfname, 'w')

    tab = '\t'

    sub_totals = {'C>A':0, 'C>G':0, 'C>T':0, 'T>A':0, 'T>C':0, 'T>G':0}

    sub_totals_list = sorted(sub_totals)

    tSNS_by_class_outfile.write ("Class" + tab + "N")
    fmsdSNS_by_class_outfile.write ("Class")
    #fsdSNS_by_class_outfile.write ("Class" + tab + "N")


    # Make a header
    for x in sub_totals_list:
        tSNS_by_class_outfile.write (tab + x)
        fmsdSNS_by_class_outfile.write (tab + x + "_mean")
        fmsdSNS_by_class_outfile.write (tab + x + "_SD")
        #fsdSNS_by_class_outfile.write (tab + x)
        fmsdSNS_by_class_outfile.write (tab + "N")

    tSNS_by_class_outfile.write ('\n')
    fmsdSNS_by_class_outfile.write ('\n')
    #fsdSNS_by_class_outfile.write ('\n')

    classes = []
    cell_classes = {}
    bio_classes = {}

    # Populate the cell_classes dict,
    # which has cells as keys and
    # classes as values
    # Use this to look up the class for
    # every cell.

    for line in class_key_file:
        if "#Sample" not in line:
            sline = line.split('\t')
            cell = sline[0]
            info_col = int(info_column)
            bio_class = sline[info_col].strip('\n')

            if bio_class not in classes:
                classes.append(bio_class)

            if bio_class not in bio_classes:
                bio_classes[bio_class] = [cell]
            else:
                bio_classes[bio_class].append(cell)

            cell_classes[cell] = bio_class

    cells_per_bio_class = {}

    for h in bio_classes:
        cells_per_bio_class[h] = str(len(bio_classes[h]))


    # Prepare the class_counter dictionary
    # which will have classes as keys and
    # a list of ints reflecting the count
    # of each SNV type as values.
    # Also prep the class_fraction dictionary
    # which will have classes as keys and

    class_counter = {}
    class_fractions = {}
    for c in classes:
        class_counter[c] = [0,0,0,0,0,0]
        class_fractions[c] = {'C>A':[], 'C>G':[], 'C>T':[], 'T>A':[], 'T>C':[], 'T>G':[]}

    # Populate the class_counter dictionary
    # For each line, if the cell is in cell_classes
    # look up it's class and store it as "SNV_class"
    # Then, in the class counter dict, at add the first
    # number in the line to the 0th positon in the
    # class counter dict.

    for line in counted_file:
        if "C>A" not in line:
            sline = line.split('\t')
            cell = sline[0]

            if cell in cell_classes:

                SNV_class = cell_classes[cell]

                for i in range (len(sline) -1):
                    class_counter[SNV_class][i] += int(sline[i+1].strip('\n'))

    for n in class_counter:
        tSNS_by_class_outfile.write(n[:-1])
        for q in class_counter[n]:
            tSNS_by_class_outfile.write(tab + str(q))

        tSNS_by_class_outfile.write('\n')

    # Populate the class_fractions dictionary
    # For each line, if the cell is in cell_classes
    # look up it's class and store it as "SNV_class"
    # Then, in the class_fractions dict, at SNV_class
    # add the first number in the line to the
    # "C>A" entry class_fractions dict,, add
    # the second number to the "C>G" entry, etc.


    for line in fraction_file:
        if "C>A" not in line:
            sline = line.split('\t')
            cell = sline[0]

            if cell in cell_classes:
                SNV_class = cell_classes[cell]
                class_fractions[SNV_class]["C>A"].append(sline[1])
                class_fractions[SNV_class]["C>G"].append(sline[2])
                class_fractions[SNV_class]["C>T"].append(sline[3])
                class_fractions[SNV_class]["T>A"].append(sline[4])
                class_fractions[SNV_class]["T>C"].append(sline[5])
                class_fractions[SNV_class]["T>G"].append(sline[6].strip('\n'))

    # Use numpy to caluclate mean and SD in the lists
    # that are currently values in the class_fractions
    # dictionary. Record the data in the class_fraction_means
    # and class_fraction_SDs dictionaries, where the keys are
    # classes, and the value for each class is a dictionary
    # with keys as SNV types and values as floats corresponding
    # to either means or SDs.

    class_fraction_means = {}
    class_fraction_SDs = {}

    for a in class_fractions:
        class_fraction_means[a] = {}
        class_fraction_SDs[a] = {}

        for snv in class_fractions[a]:
            all_snv_fractions = []
            for g in class_fractions[a][snv]:
                all_snv_fractions.append(float(g))

            class_snv_fractions = np.array(all_snv_fractions)
            class_snv_mean = np.mean(class_snv_fractions)
            class_snv_sd = np.std(class_snv_fractions)
            class_fraction_means[a][snv] = class_snv_mean
            class_fraction_SDs[a][snv] = class_snv_sd

    sorted_classes = sorted(class_fraction_means)

    for b in sorted_classes:

        fmsdSNS_by_class_outfile.write(b.strip('\r'))
        #fsdSNS_by_class_outfile.write(b[:-1])

        #fmsdSNS_by_class_outfile.write(tab + cells_per_bio_class[b])
        #fsdSNS_by_class_outfile.write(tab + cells_per_bio_class[b])

        for d in sub_totals_list:

            fmsdSNS_by_class_outfile.write(tab + str(class_fraction_means[b][d]))
            fmsdSNS_by_class_outfile.write(tab + str(class_fraction_SDs[b][d]))
            fmsdSNS_by_class_outfile.write(tab + cells_per_bio_class[b])

        fmsdSNS_by_class_outfile.write('\n')


    fraction_file.close()
    counted_file.close()
    class_key_file.close()
    tSNS_by_class_outfile.close()
    fmsdSNS_by_class_outfile.close()

def run_all():


    #fnames = glob.glob(os.path.join(dirname, "*mutations.txt"))
    # fname = "/groups/walsh/indData/Lodato/Aging/final_mutations_20170619.txt"
    print ("Running on ")
    #for fname in fnames:
    fname="/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/strandbias/ALL_mutations_AF30.bed"
    print (fname)

    #for fname in fnames:
    t_outfname = "ALL_mutations_AF30_total_SNS_by_cells"
    print (t_outfname)
    f_outfname = "ALL_mutations_AF30_fraction_SNS_by_cells"
    count_sub_types(fname, t_outfname, f_outfname)

    counted_fname="ALL_mutations_AF30_total_SNS_by_cells"
    fraction_fname="ALL_mutations_AF30_fraction_SNS_by_cells"

    class_key_fname = "/n/scratch2/eam63/merging_cells_project/merging/gvcf_files/strandbias/samples_byclass.txt"
    fmsdSNS_by_class_outfname = "%s_class_mean_SD" %(fraction_fname[:-22]) +".txt"
    tSNS_by_class_outfname = "%s_class" %(counted_fname[:-19]) +".txt"


    info_column = "1" # this info column argument is if we have several categories our class file.
    SNS_grouper(counted_fname, fraction_fname, class_key_fname, info_column,
                tSNS_by_class_outfname, fmsdSNS_by_class_outfname)




def main():
    # Main pipeline

    run_all()

if __name__ == "__main__":
    main()
