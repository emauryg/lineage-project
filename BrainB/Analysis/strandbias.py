#####################################################################################################
##
## Python code to get strand, adapted from Michael Lodato
##
####################################################################################################
import os
import glob
import numpy as np
import copy as cp


def get_strand(fname, outfname):
    """
    Get the strand information from a bedtools intersect output
    of a bed file of SNV postions and a bed file of table browser
    genes downloaded from UCSC table browser.
    """

    file = open(fname)
    outfile = open(outfname, 'w')
    strands = {}
    tab = '\t'

    # For each line, record the SNV position and the strand.
    # populate the strands dictionary with keys as SNVs and
    # lists containing strand information at index 0
    # and the entire line at index 1, as the values.
    # SNVs with more than one strand result in keys with more
    # than one value.

    for line in file:
        sline = line.split('\t')
        global_ID = sline[5]
        strand = sline[12]
        info = sline
        cell=sline[6]

        if global_ID not in strands:
            strands[global_ID] = [[strand], sline]

        else:
            strands[global_ID][0].append(strand)

    # for each SNV, write it and its strand information
    # to the outfile.  also write the other information about
    # the SNV, eg SNS type, to the outfile

    for SNV in strands:
        outfile.write(SNV)
        outfile.write('\t')
        for o in strands[SNV][0]:
            outfile.write(o)
        outfile.write('\t')
        for q in strands[SNV][1]:
            outfile.write(q.strip('\n'))
            outfile.write('\t')
        outfile.write('\n')

    file.close()
    outfile.close()

def strand_counter(strands_fname, t_outfname, f_outfname):
    """
    Use the output file of get_strand to classify SNVs
    according to strand (all + or all -), and by SNS
    type
    """

    tab = '\t'

    strands_file = open(strands_fname)
    t_outfile = open(t_outfname, "a")
    f_outfile = open(f_outfname, "a")

    cell_strands = {}
    # The cell_strands dictionary will have global_IDs as keys,
    # and for values each global ID will have
    # a list of 2 dictionaries.
    # The dictionary at index 0 will have the counts for all +
    # strand SNVs, and the dictionary at index 1
    # will ahve the counts for all - strand SNVs.

    sorted_keys = ['A>C', 'A>G', 'A>T', 'G>A', 'G>C', 'G>T']

    strand_rev_comps = {'T>G': 'A>C', 'T>C': 'A>G', 'T>A': 'A>T', 'C>T': 'G>A', 'C>G': 'G>C', 'C>A': 'G>T'}


    totals_by_cell = {}

    for line in strands_file:
        sline = line.split('\t')

        plus_strand_SNVs = {"A>C":0,"A>G":0,"A>T":0,"G>A":0,"G>C":0,"G>T":0}

        minus_strand_SNVs = {"A>C":0,"A>G":0,"A>T":0,"G>A":0,"G>C":0,"G>T":0}

        cell_ID = sline[8]
        directions = sline[1]
        SNS = sline[5] + '>' + sline[6]

        if cell_ID not in cell_strands:
            cell_strands[cell_ID] = [cp.deepcopy(plus_strand_SNVs), cp.deepcopy(minus_strand_SNVs)]
            totals_by_cell[cell_ID] = 0

        if '+' in directions:
            if '-' not in directions:
                if SNS in strand_rev_comps:
                    cell_strands[cell_ID][1][strand_rev_comps[SNS]] += 1
                else:
                    cell_strands[cell_ID][0][SNS] += 1

        elif '-' in directions:
            if '+' not in directions:
                if SNS in strand_rev_comps:
                    cell_strands[cell_ID][0][strand_rev_comps[SNS]] += 1
                else:
                    cell_strands[cell_ID][1][SNS] += 1


            totals_by_cell[cell_ID] += 1

    cell_strands_sorted_keys = sorted(cell_strands)


    for c in cell_strands_sorted_keys:
        t_outfile.write(c + "+" + tab)
        f_outfile.write(c + "+" + tab)

        for k in sorted_keys:
            t_outfile.write(str(cell_strands[c][0][k]) + tab)
            f_outfile.write(str((float(cell_strands[c][0][k]))/(float(totals_by_cell[c]))) + tab)

        t_outfile.write('\n')
        f_outfile.write('\n')

        t_outfile.write(c + "-" + tab)
        f_outfile.write(c + "-" + tab)

        for k in sorted_keys:
            t_outfile.write(str(cell_strands[c][1][k]) + tab)
            f_outfile.write(str((float(cell_strands[c][1][k]))/(float(totals_by_cell[c]))) + tab)
        t_outfile.write('\n')
        f_outfile.write('\n')

    strands_file.close()
    t_outfile.close()
    f_outfile.close()

def strand_grouper(counted_fname, fraction_fname, class_key_fname,
                tstrand_by_class_outfname, fmsdstrand_by_class_outfname):
    """
    Group strand quantified data by class, such as age group
    tissue, and disease
    """

    fraction_file = open(fraction_fname)
    counted_file = open(counted_fname)
    class_key_file = open(class_key_fname)
    tstrand_by_class_outfile = open(tstrand_by_class_outfname, 'w')
    fmsdstrand_by_class_outfile = open(fmsdstrand_by_class_outfname, 'w')

    tab = '\t'

    sub_totals = {'G>T':0, 'G>C':0, 'G>A':0, 'A>T':0, 'A>G':0, 'A>C':0}

    sub_totals_list = sorted(sub_totals)

    tstrand_by_class_outfile.write ("Class")
    fmsdstrand_by_class_outfile.write ("Class")

    # Make a header
    for x in sub_totals_list:
        tstrand_by_class_outfile.write (tab + x)

    tstrand_by_class_outfile.write ('\n')

    classes = []
    cell_classes = {}
    bio_classes = {}

    # Populate the cell_classes dict,
    # which has cells as keys and
    # classes as values
    # Use this to look up the class for
    # every cell.

    for line in class_key_file:
        if "Class" not in line:
            sline = line.split('\t')
            cell = sline[0]
            bio_class = sline[1].strip('\n').strip('\r')

            if bio_class + "+" not in classes:
                classes.append(bio_class + "+")
                classes.append(bio_class + "-")

            if bio_class + "+" not in bio_classes:
                bio_classes[bio_class + "+"] = [cell]
                bio_classes[bio_class + "-"] = [cell]

            else:
                bio_classes[bio_class + "+"].append(cell)
                bio_classes[bio_class + "-"].append(cell)

            cell_classes[cell + "+"] = bio_class + "+"
            cell_classes[cell + "-"] = bio_class + "-"

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
        class_fractions[c] = {'A>C':[], 'A>G':[], 'A>T':[], 'G>A':[], 'G>C':[], 'G>T':[]}


    for line in counted_file:
        if "C>T" not in line:
            sline = line.split('\t')
            cell_strand_data = sline[0]

            if cell_strand_data in cell_classes:
                SNV_class = cell_classes[cell_strand_data]

                for m in range(len(class_counter[SNV_class])):
                    class_counter[SNV_class][m] += int(sline[m+1].strip('\n'))

    SNS_types = ['A>C', 'A>G', 'A>T', 'G>A', 'G>C', 'G>T']

    for f in sorted(class_counter):
        tstrand_by_class_outfile.write(f)

        for n in class_counter[f]:
            tstrand_by_class_outfile.write('\t')
            tstrand_by_class_outfile.write(str(n))


        tstrand_by_class_outfile.write('\n')


    for line in fraction_file:
        if "C>A" not in line:
            sline = line.split('\t')
            cell_strand_data = sline[0]

            if cell_strand_data in cell_classes:
                SNV_class = cell_classes[cell_strand_data]
                class_fractions[SNV_class]["A>C"].append(sline[1])
                class_fractions[SNV_class]["A>G"].append(sline[2])
                class_fractions[SNV_class]["A>T"].append(sline[3])
                class_fractions[SNV_class]["G>A"].append(sline[4])
                class_fractions[SNV_class]["G>C"].append(sline[5])
                class_fractions[SNV_class]["G>T"].append(sline[6].strip('\n'))

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

    # header should have SNS types in column 1
    # col 2, 3, and 4 should be (for example) AT+ mean,
    # AT+ SD, and AT+ N, and columns 5, 6, and 7 should be
    # (for example) AT- mean, AT- SD, and AT- N.
    # use the sorted_classes list to write this line
    # then, on the next line, fill in the data for each sample

    fmsdstrand_by_class_outfile.write ("SNS type")

    for o in sorted_classes:
        fmsdstrand_by_class_outfile.write (tab + o + " mean")
        fmsdstrand_by_class_outfile.write (tab + o + " SD")
        fmsdstrand_by_class_outfile.write (tab + "N")

    fmsdstrand_by_class_outfile.write ("\n")

    for d in sub_totals_list:
        fmsdstrand_by_class_outfile.write(d)
        for l in sorted_classes:
            fmsdstrand_by_class_outfile.write(tab + str(class_fraction_means[l][d]))
            fmsdstrand_by_class_outfile.write(tab + str(class_fraction_SDs[l][d]))
            fmsdstrand_by_class_outfile.write(tab + cells_per_bio_class[l])

        fmsdstrand_by_class_outfile.write ("\n")



    fraction_file.close()
    counted_file.close()
    class_key_file.close()
    tstrand_by_class_outfile.close()
    fmsdstrand_by_class_outfile.close()


def run_all():
    'Run on all files ending with "header.bed" in dirname'

    # You can run get_strand just 1 time for each new version
    # of the mutation file. Once it produces a
    # .collapsed.txt file, you can skip the get_strand
    # function by commenting it out, and only run
    # strand_counter, where you can alter the stat filter

    # This is to run get_strand, which only needs to be run
    # once per new round of mutation files

    #fnames = glob.glob(os.path.join(dirname, "*s.bed"))
    #fname = "/groups/walsh/indData/Lodato/Aging/strandbias/final_mutations_20170619.transcripts.bed"
    fname = 'ALL_mutations_formatted_AF30.bed'
    outfname = 'ALL_mutations_formatted_AF30_stranded.bed'
    get_strand(fname, outfname)

    # This is to run strand_counter, which can be run
    # multiple times on the same collapsed file
    # That comes out of get_strand


    t_outfname = "all_strand_bias_out_total.txt"
    f_outfname = "all_strand_bias_out_fractions.txt"
    t_outfile = open(t_outfname, 'w')
    f_outfile = open(f_outfname, 'w')


    SNS_types = ['A>C', 'A>G', 'A>T', 'G>A', 'G>C', 'G>T']

    for s in SNS_types:
        t_outfile.write('\t')
        t_outfile.write(s)
        f_outfile.write('\t')
        f_outfile.write(s)

    t_outfile.write('\n')
    f_outfile.write('\n')
    t_outfile.close()
    f_outfile.close()

    #strands_fnames = glob.glob(os.path.join(dirname, "*collapsed.txt"))

    strands_fname= 'ALL_mutations_formatted_AF30_stranded.bed'
    #for strands_fname in strands_fnames:
    strand_counter(strands_fname, t_outfname, f_outfname)

    # This is to run strand_grouper, which
    # will combine cells into biological classes
    # and tally all mutations

    counted_fname = "./all_strand_bias_out_total.txt"
    fraction_fname = "./all_strand_bias_out_fractions.txt"
    tstrand_by_class_outfname = "%s_total_grouped_by_class.txt" %(counted_fname[:-4])
    fmsdstrand_by_class_outfname = "%s_meanSD_grouped_by_class.txt" %(counted_fname[:-4])

    samples = ( "1465-cortex_1-neuron_MDA_3_WGSb", "1465-cortex_1-neuron_MDA_12", "1465-cortex_1-neuron_MDA_2_WGSb", "1465-cortex_1-neuron_MDA_24", "1465-cortex_1-neuron_MDA_6_WGSb", "1465-cortex_1-neuron_MDA_18", "1465-cortex_1-neuron_MDA_39", "1465-cortex_1-neuron_MDA_47", "1465-cortex_1-neuron_MDA_51_WGSb", "1465-cortex_1-neuron_MDA_20",  "1465-cortex_1-neuron_MDA_25", "1465-cortex_1-neuron_MDA_30", "1465-cortex_1-neuron_MDA_43", "1465-cortex_1-neuron_MDA_46", "1465-cortex_1-neuron_MDA_5", "1465-cortex_1-neuron_MDA_8" )
    byclass=open("samples_byclass.txt","w")
    bioclass='neuron'
    for s in samples:
      byclass.write(s+"\t"+bioclass+'\n')
    byclass.close()

    class_key_fname = "samples_byclass.txt"

    strand_grouper(counted_fname, fraction_fname, class_key_fname,
                tstrand_by_class_outfname, fmsdstrand_by_class_outfname)

def main():

    run_all()


if __name__ == "__main__":
    main()
