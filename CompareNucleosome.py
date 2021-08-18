# author: siayouyang
# analysis of nucleosomes paired between two samples
from optparse import OptionParser
import multiprocessing
import os
from collections import Counter
import numpy
from scipy.stats import ttest_ind_from_stats

def run():
    ##OptionParser
    parser = OptionParser()
    parser.add_option("--input1", dest="input1", help="nucleosome position data 1(sorted), usually be Wt/Control")
    parser.add_option("--input2", dest="input2", help="nucleosome position data 2(sorted), usually be treatment")
    parser.add_option("-p", "--processors", dest="processors", type="int", default=5, help="Maximum processors to be run in parallel, will not exceed the maximum number of chromosomes . Default %default.")
    parser.add_option("--output1", dest="output1", help="output name for input1")
    parser.add_option("--output2", dest="output2", help="output name for input2")
    parser.add_option("--overallAbsShift", dest="overall_shift", action="store_true", help="get overall nucleosomes absalute shift statistics")
    parser.add_option("--overallSpacing", dest="overall_spacing", action="store_true", help="get overall nucleosomes spacing statistics")
    parser.add_option("--overallChange", dest="overall_change", action="store_true",help="get overall nucleosomes fold change statistics")
    parser.add_option('-F', action='store', type='int', dest='filter', default='3', help='Absolute read filter; outputs only peaks with larger read count. Default %default. ')
    parser.add_option("-g", "--gene", dest="gene", help="gene data in BED format")
    parser.add_option('--pvalue', action='store', type='float', dest='p_value', default=0.01, help='define the p-value for T-test for possible significant nucleosome shift. Default %default. ')
    parser.add_option('--occchange', action='store', type='float', dest='occchange', default=2, help='occupancy change threshold to be considered as significant,eg. 0.5 or 2 will get the same result. Default %default. ')
    parser.add_option("--significantShiftGenes", dest="sig_shift_genes", action="store_true", help="get BED with at least one significant shift nucleosome ")
    parser.add_option("--occupancyChangeGenes", dest="occ_change_genes", action="store_true", help="get BED with at least one significant occupancy change nucleosome ")

    (options, args) = parser.parse_args()


    global input1, input2, max_processors, output1, output2, overall_shift, overall_spacing, overall_change, filter, gene, p_value, occchange, sig_shift_genes, occ_change_genes
    input1 = options.input1
    input2 = options.input2
    max_processors = options.processors
    output1 = options.output1
    output2 = options.output2
    overall_shift = options.overall_shift
    overall_spacing = options.overall_spacing
    overall_change = options.overall_change
    filter = options.filter
    p_value = options.p_value
    occchange = options.occchange
    sig_shift_genes = options.sig_shift_genes
    occ_change_genes = options.occ_change_genes
    gene = options.gene

    #get chromosomes name
    chr_set = set()
    input1_file = open(input1, 'r')
    input2_file = open(input2, 'r')
    for a in input1_file:
        chr_set.add(a.split()[0])
    for b in input2_file:
        chr_set.add(b.split()[0])
    #create chromosomes name tuple
    global chr_tuple
    chr_tuple = tuple(chr_set)
    del chr_set
    input1_file.close()
    input2_file.close()

def split_chr(chr, input1, input2):
    input1_file = open(input1, 'r')
    input2_file = open(input2, 'r')
    input1_chr_temp = open(input1 + "." + chr + ".temp", 'w')
    input2_chr_temp = open(input2 + "." + chr + ".temp", 'w')
    for a in input1_file:
        if str(a.split()[0]) == str(chr) and int(a.split()[6]) >= int(filter):
            input1_chr_temp.write(a)
        else:
            pass
    for b in input2_file:
        if str(b.split()[0]) == str(chr) and int(b.split()[6]) >= int(filter):
            input2_chr_temp.write(b)
        else:
            pass
    input1_file.close()
    input2_file.close()
    input1_chr_temp.close()
    input2_chr_temp.close()


def find_pair(chr, input1, input2):
    input1_chr_temp = open(input1 + "." + chr + ".temp", 'r')
    input2_chr_temp = open(input2 + "." + chr + ".temp", 'r')
    input1_pos = []
    input2_pos = []
    for a in input1_chr_temp:
        input1_pos.append(a.split()[1])
    for b in input2_chr_temp:
        input2_pos.append(b.split()[1])
    input1_sorted_pos = sorted(input1_pos, key=int)
    input2_sorted_pos = sorted(input2_pos, key=int)
    del input1_pos, input2_pos
    input1_chr_temp.close()
    input2_chr_temp.close()
    #set
    set_1 = set()
    set_2 = set()
    if len(input1_sorted_pos) <= 1 or len(input2_sorted_pos) <= 1:
        print(f"less than 2 nucleosome in {chr} to perform comparison\n")
    else:
        for c in range(0, len(input1_sorted_pos)):
            for d in range(1, len(input2_sorted_pos)):
                shift1 = abs(int(input1_sorted_pos[c]) - int(input2_sorted_pos[d-1]))
                shift2 = abs(int(input1_sorted_pos[c]) - int(input2_sorted_pos[d]))
                if shift2 >= shift1:
                    set_1.add(f"{input1_sorted_pos[c]}:{input2_sorted_pos[d-1]}")
                    break
                else:
                    if d == len(input2_sorted_pos) - 1:
                        set_1.add(f"{input1_sorted_pos[c]}:{input2_sorted_pos[d-1]}")
                    else:
                        continue
        for c in range(0, len(input2_sorted_pos)):
            for d in range(1, len(input1_sorted_pos)):
                shift1 = abs(int(input2_sorted_pos[c]) - int(input1_sorted_pos[d-1]))
                shift2 = abs(int(input2_sorted_pos[c]) - int(input1_sorted_pos[d]))
                if shift2 >= shift1:
                    set_2.add(f"{input1_sorted_pos[d-1]}:{input2_sorted_pos[c]}")
                    break
                else:
                    if d == len(input1_sorted_pos)-1:
                        set_2.add(f"{input1_sorted_pos[d - 1]}:{input2_sorted_pos[c]}")
                    else:
                        continue
        set_1_only = set_1 - set_2
        set_2_only = set_2 - set_1
        input1_only = tuple(set_1_only)
        input2_only = tuple(set_2_only)
        input1_unpaired_pos = ()
        input2_unpaired_pos = ()
        for e in input1_only:
            input1_unpaired_pos += (int(e.split(":")[0]),)
        for f in input2_only:
            input2_unpaired_pos += (int(f.split(":")[1]),)
        del input1_sorted_pos, input2_sorted_pos, set_1, set_2, set_1_only, set_2_only, input1_only, input2_only
        #write paired chromosomes position
        input1_chr_temp = open(input1 + "." + chr + ".temp", 'r')
        input2_chr_temp = open(input2 + "." + chr + ".temp", 'r')
        input1_chr_txt = open(input1 + "." + chr + ".txt", 'w')
        input2_chr_txt = open(input2 + "." + chr + ".txt", 'w')
        for a in input1_chr_temp:
            if int(a.split()[1]) not in input1_unpaired_pos:
                input1_chr_txt.write(a)
        for b in input2_chr_temp:
            if int(b.split()[1]) not in input2_unpaired_pos:
                input2_chr_txt.write(b)
        del input1_unpaired_pos, input2_unpaired_pos
        input1_chr_temp.close()
        input2_chr_temp.close()
        input1_chr_txt.close()
        input2_chr_txt.close()
    os.remove(input1 + "." + chr + ".temp")
    os.remove(input2 + "." + chr + ".temp")


def merge_txt():
    output1_txt = open(output1, 'w')
    output2_txt = open(output2, 'w')
    for chr in chr_tuple:
        input1_chr_txt = open(input1 + "." + chr + ".txt", 'r')
        input2_chr_txt = open(input2 + "." + chr + ".txt", 'r')
        for a in input1_chr_txt:
            output1_txt.write(a)
        for b in input2_chr_txt:
            output2_txt.write(b)
        input1_chr_txt.close()
        input2_chr_txt.close()
        #os.remove(input1 + "." + chr + ".txt")
        #os.remove(input2 + "." + chr + ".txt")
    output1_txt.close()
    output2_txt.close()


def spacing():
    output1_txt = open(output1, 'r')
    output2_txt = open(output2, 'r')
    stats_txt = open(output2 + ".peaktopeak_stats.txt", 'w')
    distance_list = []
    def append_distance(output_txt):
        nucleosomes_list = []
        for a in output_txt:
            if len(nucleosomes_list) < 2:
                nucleosomes_list.append(f"{a.split()[0]}\t{a.split()[1]}")
                continue
            else:
                if str(nucleosomes_list[0].split()[0]) == str(nucleosomes_list[1].split()[0]):
                    distance = abs(int(nucleosomes_list[1].split()[1])-int(nucleosomes_list[0].split()[1]))
                    distance_list.append(distance)
                else:
                    pass
            nucleosomes_list.pop(0)
            nucleosomes_list.append(f"{a.split()[0]}\t{a.split()[1]}")
        if str(nucleosomes_list[0].split()[0]) == str(nucleosomes_list[1].split()[0]):
            distance = abs(int(nucleosomes_list[1].split()[1]) - int(nucleosomes_list[0].split()[1]))
            distance_list.append(distance)
        else:
            pass
    append_distance(output1_txt)
    append_distance(output2_txt)
    output1_txt.close()
    output2_txt.close()
    distance_count = Counter(distance_list)
    #write stats
    sample_size = int(len(distance_list))
    distance_std = numpy.std(distance_list)
    distance_mean = numpy.mean(distance_list)
    distance_min = numpy.min(distance_list)
    distance_25_percentile =numpy.quantile(distance_list, 0.25)
    distance_median = numpy.median(distance_list)
    distance_75_percentile = numpy.quantile(distance_list, 0.75)
    distance_max = numpy.max(distance_list)
    stats_txt.write(f"#peak to peak distance statistics for {output1_txt} and {output2_txt}:\n"
                    f"#sample size : {sample_size}\n"
                    f"#standard deviation :{distance_std} bp\n"
                    f"#mean : {distance_mean} bp\n"
                    f"#minimum : {distance_min} bp\n"
                    f"#25th percentile : {distance_25_percentile} bp\n"
                    f"#median : {distance_median} bp\n"
                    f"#75th percentile : {distance_75_percentile} bp\n"
                    f"#maximum : {distance_max} bp\n"
                    f"#------------------------------------\n"
                    f"#distance\tcount\n")
    for b in sorted(distance_count, key = int):
        stats_txt.write(f"{int(b)}\t{int(distance_count[b])}\n")
    stats_txt.close()
    del distance_list, distance_count


def absolute_shift():
    output1_txt = open(output1, 'r')
    output2_txt = open(output2, 'r')
    stats_txt = open(output2 + ".absshift_stats.txt", 'w')
    abs_shift_list = []
    output1_pos = []
    output2_pos = []
    for a in output1_txt:
        output1_pos.append(a.split()[1])
    for b in output2_txt:
        output2_pos.append(b.split()[1])
    for c in range(0, len(output1_pos)):
        abs_shift = abs(int(output1_pos[c])-int(output2_pos[c]))
        abs_shift_list.append(abs_shift)
    output1_txt.close()
    output2_txt.close()
    del output1_pos, output2_pos
    shift_count = Counter(abs_shift_list)
    # write stats
    sample_size = int(len(abs_shift_list))
    shift_std = numpy.std(abs_shift_list)
    shift_mean = numpy.mean(abs_shift_list)
    shift_min = numpy.min(abs_shift_list)
    shift_25_percentile = numpy.quantile(abs_shift_list, 0.25)
    shift_median = numpy.median(abs_shift_list)
    shift_75_percentile = numpy.quantile(abs_shift_list, 0.75)
    shift_max = numpy.max(abs_shift_list)
    stats_txt.write(f"#absolute nucleosome shift statistics between {output1_txt} and {output2_txt}:\n"
                    f"#sample size : {sample_size}\n"
                    f"#standard deviation :{shift_std} bp\n"
                    f"#mean : {shift_mean} bp\n"
                    f"#minimum : {shift_min} bp\n"
                    f"#25th percentile : {shift_25_percentile} bp\n"
                    f"#median : {shift_median} bp\n"
                    f"#75th percentile : {shift_75_percentile} bp\n"
                    f"#maximum : {shift_max} bp\n"
                    f"#------------------------------------\n"
                    f"#absolute shift\tcount\n")
    for d in sorted(shift_count, key = int):
        stats_txt.write(f"{int(d)}\t{int(shift_count[d])}\n")
    stats_txt.close()
    del abs_shift_list, shift_count


def occupancy_change():
    output1_txt = open(output1, 'r')
    output2_txt = open(output2, 'r')
    stats_txt = open(output2 + ".occchange_stats.txt", 'w')
    occ_change_list = []
    output1_occ = []
    output2_occ = []
    for a in output1_txt:
        output1_occ.append(a.split()[6])
    for b in output2_txt:
        output2_occ.append(b.split()[6])
    for c in range(0, len(output1_occ)):
        occ_change = float(int(output2_occ[c])/int(output1_occ[c]))
        occ_change_list.append(occ_change)
    output1_txt.close()
    output2_txt.close()
    del output1_occ, output2_occ
    occ_change_count = Counter(occ_change_list)
    # write stats
    sample_size = int(len(occ_change_list))
    occ_change_std = numpy.std(occ_change_list)
    occ_change_mean = numpy.mean(occ_change_list)
    occ_change_min = numpy.min(occ_change_list)
    occ_change_25_percentile = numpy.quantile(occ_change_list, 0.25)
    occ_change_median = numpy.median(occ_change_list)
    occ_change_75_percentile = numpy.quantile(occ_change_list, 0.75)
    occ_change_max = numpy.max(occ_change_list)
    stats_txt.write(f"#nucleosome occupancy change statistics between {output2_txt}/{output1_txt}:\n"
                    f"#sample size : {sample_size}\n"
                    f"#standard deviation :{occ_change_std}\n"
                    f"#mean : {occ_change_mean}\n"
                    f"#minimum : {occ_change_min}\n"
                    f"#25th percentile : {occ_change_25_percentile}\n"
                    f"#median : {occ_change_median}\n"
                    f"#75th percentile : {occ_change_75_percentile}\n"
                    f"#maximum : {occ_change_max}\n"
                    f"#------------------------------------\n"
                    f"#fold change\tcount\n")
    for d in sorted(occ_change_count, key=float):
        stats_txt.write(f"{float(d)}\t{int(occ_change_count[d])}\n")
    stats_txt.close()
    del occ_change_list, occ_change_count


def significant_shift_genes(chr, input1, input2):
    # genes >= 560bp
    def filter_gene_length(chr):
        gene_raw_file = open(gene, 'r')
        gene_chr_file = open(gene + chr + ".clean.txt", 'w')
        for a in gene_raw_file:
            if (abs(int(a.split()[2]) - int(a.split()[1])) >= 560) and str(a.split()[0]) == str(chr):
                gene_chr_file.write(a)
            else:
                pass
        gene_raw_file.close()
        gene_chr_file.close()
    filter_gene_length(chr)

    def shift_pvalue(chr, input1, input2):
        input1_chr_txt = open(input1 + "." + chr + ".txt", 'r')
        input2_chr_txt = open(input2 + "." + chr + ".txt", 'r')
        input1_pos_tuple = ()
        input2_pos_tuple = ()
        input1_std_tuple = ()
        input2_std_tuple = ()
        input1_mean_tuple = ()
        input2_mean_tuple = ()
        input1_readcount_tuple = ()
        input2_readcount_tuple = ()
        for a in input1_chr_txt:
            input1_pos_tuple += (int(a.split()[1]),)
            input1_std_tuple += (float(a.split()[4]),)
            input1_mean_tuple += (float(a.split()[5]),)
            input1_readcount_tuple += (int(a.split()[6]),)
        for b in input2_chr_txt:
            input2_pos_tuple += (int(b.split()[1]),)
            input2_std_tuple += (float(b.split()[4]),)
            input2_mean_tuple += (float(b.split()[5]),)
            input2_readcount_tuple += (int(b.split()[6]),)
        input1_chr_txt.close()
        input2_chr_txt.close()
        input1_chr_ttest_temp = open(output1 + "." + chr + ".ttest.temp", 'w')
        for c in range(0, len(input1_pos_tuple)):
            pos = int(input1_pos_tuple[c])
            shift = int(input2_pos_tuple[c] - input1_pos_tuple[c])
            mean1 = float(input1_mean_tuple[c])
            mean2 = float(input2_mean_tuple[c])
            std1 = float(input1_std_tuple[c])
            std2 = float(input2_std_tuple[c])
            readcount1 = int(input1_readcount_tuple[c])
            readcount2 = int(input2_readcount_tuple[c])
            ttest = ttest_ind_from_stats(mean1, std1, readcount1, mean2, std2, readcount2, equal_var=False, alternative='two-sided')
            pvalue = float(ttest[1])
            occ_change = float(readcount2/readcount1)
            input1_chr_ttest_temp.write(f"{chr}\t{pos}\t{shift}\t{pvalue}\t{occ_change}\n")
        input1_chr_ttest_temp.close()
        del input1_pos_tuple, input2_pos_tuple, input1_std_tuple, input2_std_tuple, input1_mean_tuple, input2_mean_tuple, input1_readcount_tuple, input2_readcount_tuple
    shift_pvalue(chr, input1, input2)

    def relate_to_gene(chr):
        gene_chr_file = open(gene + chr + ".clean.txt", 'r')
        input1_chr_ttest_temp = open(output1 + "." + chr + ".ttest.temp", 'r')
        input1_chr_ttest_readlines = input1_chr_ttest_temp.readlines()
        tss4nuc_chr_file = open(output1 + chr + ".TSS4nuc.sigshift.bed", 'w')
        if len(input1_chr_ttest_readlines) <= 1:
            print(f"less than 2 nucleosome in {chr} to perform comparison\n")
        for a in gene_chr_file:
            for b in range(1, len(input1_chr_ttest_readlines)):
                if str(a.split()[5]) == "+":
                    previous_distance = abs(int(a.split()[1]) - int(input1_chr_ttest_readlines[b-1].split()[1]))
                    current_distance = abs(int(a.split()[1]) - int(input1_chr_ttest_readlines[b].split()[1]))
                    if (current_distance >= previous_distance) and (previous_distance <= 150):
                        if (b - 1) >= 0 and ((b - 1) < int(len(input1_chr_ttest_readlines))):
                            plus1_pos = int(input1_chr_ttest_readlines[b - 1].split()[1])
                            plus1_shift = int(input1_chr_ttest_readlines[b-1].split()[2])
                            plus1_pvalue = float(input1_chr_ttest_readlines[b - 1].split()[3])
                            plus1_occ_change = float(input1_chr_ttest_readlines[b - 1].split()[4])
                            if (b) >= 0 and ((b) < int(len(input1_chr_ttest_readlines))):
                                plus2_pos = int(input1_chr_ttest_readlines[b].split()[1])
                                plus2_shift = int(input1_chr_ttest_readlines[b].split()[2])
                                plus2_pvalue = float(input1_chr_ttest_readlines[b].split()[3])
                                plus2_occ_change = float(input1_chr_ttest_readlines[b].split()[4])
                                if (b + 1) >= 0 and ((b + 1) < int(len(input1_chr_ttest_readlines))):
                                    plus3_pos = int(input1_chr_ttest_readlines[b+1].split()[1])
                                    plus3_shift = int(input1_chr_ttest_readlines[b+1].split()[2])
                                    plus3_pvalue = float(input1_chr_ttest_readlines[b + 1].split()[3])
                                    plus3_occ_change = float(input1_chr_ttest_readlines[b + 1].split()[4])
                                    if (b + 2) >= 0 and ((b + 2) < int(len(input1_chr_ttest_readlines))):
                                        plus4_pos = int(input1_chr_ttest_readlines[b + 2].split()[1])
                                        plus4_shift = int(input1_chr_ttest_readlines[b + 2].split()[2])
                                        plus4_pvalue = float(input1_chr_ttest_readlines[b + 2].split()[3])
                                        plus4_occ_change = float(input1_chr_ttest_readlines[b + 2].split()[4])
                                        if abs(plus1_pos - plus2_pos) <= 250 and abs(plus2_pos - plus3_pos) <= 250 and abs(plus3_pos - plus4_pos) <= 250:
                                            tss4nuc_chr_file.write(f"{chr}\t{plus1_pos}\t{plus1_pos + 1}\t{a.split()[3]}\t{0}\t{a.split()[5]}\t{plus1_shift}\t{plus1_pvalue}\t{plus2_shift}\t{plus2_pvalue}\t{plus3_shift}\t{plus3_pvalue}\t{plus4_shift}\t{plus4_pvalue}\t{plus1_occ_change}\t{plus2_occ_change}\t{plus3_occ_change}\t{plus4_occ_change}\n")
                                            break
                                        else:
                                            break
                                    else:
                                        break
                                else:
                                    break
                            else:
                                break
                        else:
                            break
                elif str(a.split()[5]) == "-":
                    previous_distance = abs(int(a.split()[2]) - int(input1_chr_ttest_readlines[b - 1].split()[1]))
                    current_distance = abs(int(a.split()[2]) - int(input1_chr_ttest_readlines[b].split()[1]))
                    if (current_distance >= previous_distance) and (previous_distance <= 150):
                        if (b - 1) >= 0 and ((b - 1) < int(len(input1_chr_ttest_readlines))):
                            plus1_pos = int(input1_chr_ttest_readlines[b - 1].split()[1])
                            plus1_shift = -1 * int(input1_chr_ttest_readlines[b - 1].split()[2])
                            plus1_pvalue = float(input1_chr_ttest_readlines[b - 1].split()[3])
                            plus1_occ_change = float(input1_chr_ttest_readlines[b - 1].split()[4])
                            if (b - 2) >= 0 and ((b - 2) < int(len(input1_chr_ttest_readlines))):
                                plus2_pos = int(input1_chr_ttest_readlines[b - 2].split()[1])
                                plus2_shift = -1 * int(input1_chr_ttest_readlines[b - 2].split()[2])
                                plus2_pvalue = float(input1_chr_ttest_readlines[b - 2].split()[3])
                                plus2_occ_change = float(input1_chr_ttest_readlines[b - 2].split()[4])
                                if (b - 3) >= 0 and ((b - 3) < int(len(input1_chr_ttest_readlines))):
                                    plus3_pos = int(input1_chr_ttest_readlines[b - 3].split()[1])
                                    plus3_shift = -1 * int(input1_chr_ttest_readlines[b - 3].split()[2])
                                    plus3_pvalue = float(input1_chr_ttest_readlines[b - 3].split()[3])
                                    plus3_occ_change = float(input1_chr_ttest_readlines[b - 3].split()[4])
                                    if (b - 4) >= 0 and ((b - 4) < int(len(input1_chr_ttest_readlines))):
                                        plus4_pos = int(input1_chr_ttest_readlines[b - 4].split()[1])
                                        plus4_shift = -1 * int(input1_chr_ttest_readlines[b - 4].split()[2])
                                        plus4_pvalue = float(input1_chr_ttest_readlines[b - 4].split()[3])
                                        plus4_occ_change = float(input1_chr_ttest_readlines[b - 4].split()[4])
                                        if abs(plus1_pos - plus2_pos) <= 250 and abs(plus2_pos - plus3_pos) <= 250 and abs(plus3_pos - plus4_pos) <= 250:
                                            tss4nuc_chr_file.write(f"{chr}\t{plus1_pos - 1}\t{plus1_pos}\t{a.split()[3]}\t{0}\t{a.split()[5]}\t{plus1_shift}\t{plus1_pvalue}\t{plus2_shift}\t{plus2_pvalue}\t{plus3_shift}\t{plus3_pvalue}\t{plus4_shift}\t{plus4_pvalue}\t{plus1_occ_change}\t{plus2_occ_change}\t{plus3_occ_change}\t{plus4_occ_change}\n")
                                            break
                                        else:
                                            break
                                    else:
                                        break
                                else:
                                    break
                            else:
                                break
                        else:
                            break
                else:
                    print("gene data[5] format error")
                    break
        gene_chr_file.close()
        input1_chr_ttest_temp.close()
        tss4nuc_chr_file.close()
        del input1_chr_ttest_readlines
    relate_to_gene(chr)
    os.remove(output1 + "." + chr + ".ttest.temp")
    os.remove(gene + chr + ".clean.txt")


def merge_sig_shift_genes():
    tss4nuc_file = open(output1 + ".TSS4nuc.sigshift.bed", 'w')
    for chr in chr_tuple:
        tss4nuc_chr_file = open(output1 + chr + ".TSS4nuc.sigshift.bed", 'r')
        for a in tss4nuc_chr_file:
            plus1_pvalue = float(a.split()[7])
            plus2_pvalue = float(a.split()[9])
            plus3_pvalue = float(a.split()[11])
            plus4_pvalue = float(a.split()[13])
            if plus1_pvalue <= p_value or plus2_pvalue <= p_value or plus3_pvalue <= p_value or plus4_pvalue <= p_value:
                tss4nuc_file.write(a)
            else:
                pass
        tss4nuc_chr_file.close()
    tss4nuc_file.close()
    # get info
    tss4nuc_file = open(output1 + ".TSS4nuc.sigshift.bed", 'r')
    tss4nuc_info_file = open(output1 + ".TSS4nuc.sigshift.info.txt", 'w')
    plus1_shift_list = []
    plus2_shift_list = []
    plus3_shift_list = []
    plus4_shift_list = []
    for a in tss4nuc_file:
        plus1_shift_list.append(int(a.split()[6]))
        plus2_shift_list.append(int(a.split()[8]))
        plus3_shift_list.append(int(a.split()[10]))
        plus4_shift_list.append(int(a.split()[12]))
    plus1_mean = numpy.mean(plus1_shift_list)
    plus1_median = numpy.median(plus1_shift_list)
    plus1_min = numpy.min(plus1_shift_list)
    plus1_max = numpy.max(plus1_shift_list)
    plus1_std = numpy.std(plus1_shift_list)
    plus1_shift_count =  Counter(plus1_shift_list)
    plus2_mean = numpy.mean(plus2_shift_list)
    plus2_median = numpy.median(plus2_shift_list)
    plus2_min = numpy.min(plus2_shift_list)
    plus2_max = numpy.max(plus2_shift_list)
    plus2_std = numpy.std(plus2_shift_list)
    plus2_shift_count = Counter(plus2_shift_list)
    plus3_mean = numpy.mean(plus3_shift_list)
    plus3_median = numpy.median(plus3_shift_list)
    plus3_min = numpy.min(plus3_shift_list)
    plus3_max = numpy.max(plus3_shift_list)
    plus3_std = numpy.std(plus3_shift_list)
    plus3_shift_count = Counter(plus3_shift_list)
    plus4_mean = numpy.mean(plus4_shift_list)
    plus4_median = numpy.median(plus4_shift_list)
    plus4_min = numpy.min(plus4_shift_list)
    plus4_max = numpy.max(plus4_shift_list)
    plus4_std = numpy.std(plus4_shift_list)
    plus4_shift_count = Counter(plus4_shift_list)
    tss4nuc_info_file.write(f"# +1 to +4 nucleosome shift statistics:\n"
                            f"# +1 mean = {plus1_mean}\t+2 mean = {plus2_mean}\t+3 mean = {plus3_mean}\t+4 mean = {plus4_mean}\n"
                            f"# +1 median = {plus1_median}\t+2 median = {plus2_median}\t+3 median = {plus3_median}\t+4 median = {plus4_median}\n"
                            f"# +1 min = {plus1_min}\t+2 min = {plus2_min}\t+3 min = {plus3_min}\t+4 min = {plus4_min}\n"
                            f"# +1 max = {plus1_max}\t+2 max = {plus2_max}\t+3 max = {plus3_max}\t+4 max = {plus4_max}\n"
                            f"# +1 std = {plus1_std}\t+2 std = {plus2_std}\t+3 std = {plus3_std}\t+4 std = {plus4_std}\n"
                            f"# ---------------------------------------------------\n")
    max_count_len = max([int(len(plus1_shift_count)), int(len(plus2_shift_count)), int(len(plus3_shift_count)), int(len(plus4_shift_count))])
    for b in range(0, max_count_len):
        def shift_count(plus_shift_count, idx):
            if idx < len(plus_shift_count):
                shift = (sorted(plus_shift_count, key=int))[idx]
                count = plus_shift_count[shift]
            else:
                shift = "NA"
                count = "NA"
            return shift, count
        plus1_shift, plus1_count = shift_count(plus1_shift_count, b)
        plus2_shift, plus2_count = shift_count(plus2_shift_count, b)
        plus3_shift, plus3_count = shift_count(plus3_shift_count, b)
        plus4_shift, plus4_count = shift_count(plus4_shift_count, b)
        tss4nuc_info_file.write(f"{plus1_shift}\t{plus1_count}\t{plus2_shift}\t{plus2_count}\t{plus3_shift}\t{plus3_count}\t{plus4_shift}\t{plus4_count}\n")
    tss4nuc_file.close()
    tss4nuc_info_file.close()


def merge_occ_change_genes():
    tss4nuc_file = open(output1 + ".TSS4nuc.occchange.bed", 'w')
    for chr in chr_tuple:
        tss4nuc_chr_file = open(output1 + chr + ".TSS4nuc.sigshift.bed", 'r')
        for a in tss4nuc_chr_file:
            plus1_occhange = float(a.split()[14])
            plus2_occhange = float(a.split()[15])
            plus3_occhange = float(a.split()[16])
            plus4_occhange = float(a.split()[17])
            greater_occchange = float(max(float(occchange), float(1/occchange)))
            lesser_occchange = float(min(float(occchange), float(1 / occchange)))
            if ((plus1_occhange <= lesser_occchange) or (plus1_occhange >= greater_occchange)) \
                    or ((plus2_occhange <= lesser_occchange) or (plus2_occhange >= greater_occchange)) \
                    or ((plus3_occhange <= lesser_occchange) or (plus3_occhange >= greater_occchange)) \
                    or ((plus4_occhange <= lesser_occchange) or (plus4_occhange >= greater_occchange)):
                tss4nuc_file.write(a)
            else:
                pass
        tss4nuc_chr_file.close()
    tss4nuc_file.close()
    # get info
    tss4nuc_file = open(output1 + ".TSS4nuc.occchange.bed", 'r')
    tss4nuc_info_file = open(output1 + ".TSS4nuc.occchange.info.txt", 'w')
    plus1_change_list = []
    plus2_change_list = []
    plus3_change_list = []
    plus4_change_list = []
    for a in tss4nuc_file:
        plus1_change_list.append(float(a.split()[14]))
        plus2_change_list.append(float(a.split()[15]))
        plus3_change_list.append(float(a.split()[16]))
        plus4_change_list.append(float(a.split()[17]))
    plus1_mean = numpy.mean(plus1_change_list)
    plus1_median = numpy.median(plus1_change_list)
    plus1_min = numpy.min(plus1_change_list)
    plus1_max = numpy.max(plus1_change_list)
    plus1_std = numpy.std(plus1_change_list)
    plus1_change_count = Counter(plus1_change_list)
    plus2_mean = numpy.mean(plus2_change_list)
    plus2_median = numpy.median(plus2_change_list)
    plus2_min = numpy.min(plus2_change_list)
    plus2_max = numpy.max(plus2_change_list)
    plus2_std = numpy.std(plus2_change_list)
    plus2_change_count = Counter(plus2_change_list)
    plus3_mean = numpy.mean(plus3_change_list)
    plus3_median = numpy.median(plus3_change_list)
    plus3_min = numpy.min(plus3_change_list)
    plus3_max = numpy.max(plus3_change_list)
    plus3_std = numpy.std(plus3_change_list)
    plus3_change_count = Counter(plus3_change_list)
    plus4_mean = numpy.mean(plus4_change_list)
    plus4_median = numpy.median(plus4_change_list)
    plus4_min = numpy.min(plus4_change_list)
    plus4_max = numpy.max(plus4_change_list)
    plus4_std = numpy.std(plus4_change_list)
    plus4_change_count = Counter(plus4_change_list)
    tss4nuc_info_file.write(f"# +1 to +4 nucleosome occupancy change statistics:\n"
                            f"# +1 mean = {plus1_mean}\t+2 mean = {plus2_mean}\t+3 mean = {plus3_mean}\t+4 mean = {plus4_mean}\n"
                            f"# +1 median = {plus1_median}\t+2 median = {plus2_median}\t+3 median = {plus3_median}\t+4 median = {plus4_median}\n"
                            f"# +1 min = {plus1_min}\t+2 min = {plus2_min}\t+3 min = {plus3_min}\t+4 min = {plus4_min}\n"
                            f"# +1 max = {plus1_max}\t+2 max = {plus2_max}\t+3 max = {plus3_max}\t+4 max = {plus4_max}\n"
                            f"# +1 std = {plus1_std}\t+2 std = {plus2_std}\t+3 std = {plus3_std}\t+4 std = {plus4_std}\n"
                            f"# ---------------------------------------------------\n")
    max_count_len = max([int(len(plus1_change_count)), int(len(plus2_change_count)), int(len(plus3_change_count)),
                         int(len(plus4_change_count))])
    for b in range(0, max_count_len):
        def change_count(plus_change_count, idx):
            if idx < len(plus_change_count):
                change = (sorted(plus_change_count, key=float))[idx]
                count = plus_change_count[change]
            else:
                change = "NA"
                count = "NA"
            return change, count

        plus1_change, plus1_count = change_count(plus1_change_count, b)
        plus2_change, plus2_count = change_count(plus2_change_count, b)
        plus3_change, plus3_count = change_count(plus3_change_count, b)
        plus4_change, plus4_count = change_count(plus4_change_count, b)
        tss4nuc_info_file.write(
            f"{plus1_change}\t{plus1_count}\t{plus2_change}\t{plus2_count}\t{plus3_change}\t{plus3_count}\t{plus4_change}\t{plus4_count}\n")
    tss4nuc_file.close()
    tss4nuc_info_file.close()


if __name__ == "__main__":
    run()

    def multi_processes(func):
        for i in range(0, len(chr_tuple), max_processors):
            if int(i + max_processors) > int(len(chr_tuple)):
                chr_tuple_chunk = chr_tuple[i:int(len(chr_tuple))]
            else:
                chr_tuple_chunk = chr_tuple[i:i + max_processors]
            processes = []
            for chr in chr_tuple_chunk:
                processes.append(multiprocessing.Process(target=func, args=(chr,input1, input2)))
            for process in processes:
                process.start()
            for process in processes:
                process.join()


    multi_processes(split_chr)
    multi_processes(find_pair)
    merge_txt()
    if overall_spacing:
        spacing()
    if overall_shift:
        absolute_shift()
    if overall_change:
        occupancy_change()

    if sig_shift_genes or occ_change_genes:
        multi_processes(significant_shift_genes)
        if sig_shift_genes:
            merge_sig_shift_genes()
        if occ_change_genes:
            merge_occ_change_genes()
        for chr in chr_tuple:
            os.remove(output1 + chr + ".TSS4nuc.sigshift.bed")

    for chr in chr_tuple:
        os.remove(input1 + "." + chr + ".txt")
        os.remove(input2 + "." + chr + ".txt")

    print("Done!")