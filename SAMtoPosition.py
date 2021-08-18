# author: siayouyang
# for sacCer3 genome
# for center 1bp read mapping
from optparse import OptionParser
import os
import numpy
from collections import Counter
import time
import math
import multiprocessing


def run():
    ##OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="sam", help="properly paired input data in SAM format")
    parser.add_option("-s","--sigma", dest="sigma", type="float", default=20,help="Sigma to be used when smoothing. Default %default")
    parser.add_option("-e", "--exclusion", dest="exclusion", type="int", default=147, help="Exclusion zone around each peak that prevents others from being called. Default %default.")
    parser.add_option("--rawWig", dest="rawwig", action="store_true", help="save raw midpoint position in WIG format")
    parser.add_option("--smoothWig", dest="smoothwig", action="store_true", help="save smooth data in WIG format")
    parser.add_option("-p", "--processors", dest="processors", type="int", default=5, help="Maximum processors to be run in parallel, will not exceed the maximum number of chromosomes. Default %default.")
    parser.add_option('-F', action='store', type='int', dest='filter', default='3',help='Absolute read filter; outputs only peaks with larger read count. Default %default. ')
    (options, args) = parser.parse_args()

    global sam, sigma, exclusion, rawwig, smoothwig, max_processors, filter
    sam = options.sam
    sigma = options.sigma
    exclusion = options.exclusion
    rawwig = options.rawwig
    smoothwig = options.smoothwig
    max_processors = options.processors
    filter = options.filter

#chromosomes info for SacCer3
chr_name = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII",
           "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrM")
chr_len = (230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177,
               924431, 784333, 1091291, 948066, 85779)


def find_center():
    chr_set = set()
    sam_file = open(sam, 'r')
    temp_file = open(sam + ".temp", 'w')
    for a in sam_file:
        if a.startswith("@"):
            continue
        else:
            chr_set.add(a.split()[2])
            chr = a.split()[2]
            pos = int(a.split()[3])
            length = int(a.split()[8])
            if length > 0:
                temp_file.write(f'{chr}\t{pos + (length // 2) - 1}\n')
            else:
                continue
    #will only proceed chrmosomes in the tuple in the rest of the scripts
    global chr_tuple
    chr_tuple = tuple(chr_set)
    del chr_set
    sam_file.close()
    temp_file.close()


def split_chr(a, sam, chr_name, chr_len, sigma, exclusion):
    #create .temp for each chromosome
    temp_file = open(sam + ".temp", 'r')
    chr_temp_file = open(sam + "." + a + ".temp", 'w')
    for b in temp_file:
        if str(b.split()[0]) == str(a):
            chr_temp_file.write(b)
    chr_temp_file.close()
    temp_file.close()


def count_reads(a, sam, chr_name, chr_len, sigma, exclusion):
    chr_temp_file = open(sam + "." + a + ".temp", 'r')
    chr_temp2_file = open(sam + "." + a + ".2.temp", 'w')
    pos_tuple = ()
    for b in chr_temp_file:
        pos_tuple += (b.split()[1],)
    count_dict = Counter(pos_tuple)
    for c in sorted(count_dict, key=int):
        chr = a
        start = int(c)
        chr_temp2_file.write(f'{chr}\t{start}\t{count_dict[c]}\n')
    chr_temp_file.close()
    chr_temp2_file.close()
    os.remove(sam + "." + a + ".temp")
    del pos_tuple, count_dict


def create_wig(a, sam, chr_name, chr_len, sigma, exclusion):
    chr_temp2_file = open(sam + "." + a + ".2.temp", 'r')
    chr_wig_file = open(sam + "." + a + ".wig", 'w')
    start = 0
    for b in chr_temp2_file:
        end = int(b.split()[1])
        for c in range(start, end):
            chr_wig_file.write(f'{a}\t{int(c)}\t{int(c)+1}\t0\n')
        chr_wig_file.write(f'{a}\t{int(b.split()[1])}\t{int(b.split()[1]) + 1}\t{int(b.split()[2])}\n')
        start = end + 1
    for d in range(start, chr_len[chr_name.index(a)]):
            chr_wig_file.write(f'{a}\t{int(d)}\t{int(d)+1}\t0\n')
    chr_temp2_file.close()
    chr_wig_file.close()
    os.remove(sam + "." + a + ".2.temp")


def merge_raw_wig():
    raw_wig_file = open(sam + ".wig", 'w')
    for a in chr_tuple:
        chr_wig_file = open(sam + "." + a + ".wig", 'r')
        raw_wig_file = open(sam + ".wig", 'a')
        for b in chr_wig_file:
            raw_wig_file.write(b)
        chr_wig_file.close()
    raw_wig_file.close()


def gaussian_smoothing(a, sam, chr_name, chr_len, sigma, exclusion):
    width=(4 * sigma)
    def normal_func(x):
        return math.exp(-x * x / (2 * (sigma ** 2)))
    gaussian = list(map(normal_func, range(-width, width)))
    gaussian = numpy.array(gaussian, numpy.float)
    # normalization
    gaussian_list = 1.0 / math.sqrt(2 * numpy.pi * (sigma ** 2)) * gaussian
    chr_wig_file = open(sam + "." + a + ".wig", 'r')
    chr_wig_file_readlines = tuple(chr_wig_file.readlines())
    chr_smooth_file = open(sam + "." + a + ".smooth.wig", 'w')
    for b in range(0, len(chr_wig_file_readlines)):
        raw_list = []
        for c in range((int(b) - width), (int(b) + width)):
            if c < 0:
                raw_list.append(float(0))
            elif c >= (int(chr_len[chr_name.index(a)])):
                raw_list.append(float(0))
            else:
                raw_list.append(float(chr_wig_file_readlines[c].split()[3]))
        smooth_list = []
        for d in range(0, len(gaussian_list)):
            e = float(gaussian_list[d]) * float(raw_list[d])
            smooth_list.append(e)
        smooth_value = sum(smooth_list)
        chr_smooth_file.write(f'{chr_wig_file_readlines[b].split()[0]}\t{chr_wig_file_readlines[b].split()[1]}\t{chr_wig_file_readlines[b].split()[2]}\t{smooth_value}\n')
    chr_wig_file.close()
    chr_smooth_file.close()
    del raw_list, smooth_list, chr_wig_file_readlines
    del gaussian_list


def find_peaks(a, sam, chr_name, chr_len, sigma, exclusion):
    smooth_file = open(sam + "." + a + ".smooth.wig", 'r')
    peaks_file = open(sam + "." + a + ".peaks.temp", 'w')
    peak_list = []
    for b in smooth_file:
        if len(peak_list) < 3:
            peak_list.append(b.split()[3])
            continue
        else:
            if (float(peak_list[1]) > float(peak_list[0])) and (float(peak_list[1]) > float(peak_list[2])):
                val = float(peak_list[1])
                pos = int(b.split()[1]) - 2
                peaks_file.write(f'{a}\t{pos}\t{val}\n')
            else:
                pass
        peak_list.pop(0)
        peak_list.append(b.split()[3])
    if (float(peak_list[1]) > float(peak_list[0])) and (float(peak_list[1]) > float(peak_list[2])):
        val = float(peak_list[1])
        pos = int(b.split()[1]) - 1
        peaks_file.write(f'{a}\t{pos}\t{val}\n')
    smooth_file.close()
    peaks_file.close()
    del peak_list


def merge_smooth_wig():
    smooth_wig_file = open(sam + ".smooth.wig", 'w')
    for a in chr_tuple:
        smooth_file = open(sam + "." + a + ".smooth.wig", 'r')
        smooth_wig_file = open(sam + ".smooth.wig", 'a')
        for b in smooth_file:
            smooth_wig_file.write(b)
        smooth_file.close()
    smooth_wig_file.close()


def perform_exclusion(a, sam, chr_name, chr_len, sigma, exclusion):
    peaks_file = open(sam + "." + a + ".peaks.temp", 'r')
    exclude_file = open(sam + "." + a + ".exclude.temp", 'w')

    #create dictionary
    list_pos = []
    list_val = []
    for b in peaks_file:
        list_pos.append(int(b.split()[1]))
        list_val.append(float(b.split()[2]))   #peak_height
    peaks_file.close()

    dict = {}
    for c in range(0, len(list_pos)):
        dict[list_pos[c]] = list_val[c]

    safe_keys = []
    while len(dict.keys()) > 0:
        max_value = (max(dict.values()))

        def get_key(dict, value):
            return [k for k, v in dict.items() if v == value]

        max_key = get_key(dict, max_value)
        safe_keys.append(max_key[0])

        exclusion_range = numpy.arange(int(max_key[0]) - int(exclusion//2), int(max_key[0]) + int(exclusion//2))
        for r in exclusion_range:
            if r in dict.keys():
                del dict[r]
            else:
                pass
    peaks_file = open(sam + "." + a + ".peaks.temp", 'r')
    for d in peaks_file:
        if int(d.split()[1]) in safe_keys:
            exclude_file.write(d)
        else:
            pass
    peaks_file.close()
    exclude_file.close()
    os.remove(sam + "." + a + ".peaks.temp")
    os.remove(sam + "." + a + ".smooth.wig")


#std,mean,readcount
def get_info(a, sam, chr_name, chr_len, sigma, exclusion):
    exclude_file = open(sam + "." + a + ".exclude.temp", 'r')
    exclude_file_readlines = exclude_file.readlines()
    chr_wig_file = open(sam + "." + a + ".wig", 'r')
    chr_wig_tuple = ()
    for chr_wig_line in chr_wig_file:
        chr_wig_tuple += (float(chr_wig_line.split()[3]),)
    info_file = open(sam + "." + a + ".info.temp", 'w')
    for b in range(0, len(exclude_file_readlines)):
        # width=146bp
        peak_range = numpy.arange((int(exclude_file_readlines[b].split()[1]) - 73),
                                  (int(exclude_file_readlines[b].split()[1]) + 73))
        position_list = []
        frequency_list = []
        for c in peak_range:
            if c < 0:
                pass
            elif c >= (len(chr_wig_tuple)):
                pass
            else:
                position_list.append(int(c))
                frequency_list.append(float(chr_wig_tuple[c]))
        flat_list = []
        for d in range(0, len(position_list)):
            position = int(position_list[d])
            frequency = int(float(frequency_list[d]))
            n = 0
            while n < frequency:
                flat_list.append(position)
                n += 1
        read_count = int(len(flat_list))
        peaks_std = numpy.std(flat_list)
        mean = numpy.mean(flat_list)
        if read_count >= int(filter):
            info_file.write(f'{exclude_file_readlines[b].split()[0]}\t{int(exclude_file_readlines[b].split()[1])}\t{int(exclude_file_readlines[b].split()[1])+1}\t{exclude_file_readlines[b].split()[2]}\t{peaks_std}\t{mean}\t{read_count}\n')
    exclude_file.close()
    chr_wig_file.close()
    info_file.close()
    os.remove(sam + "." + a + ".exclude.temp")
    os.remove(sam + "." + a + ".wig")
    del exclude_file_readlines, chr_wig_tuple


def merge_info():
    merged_info_file = open(sam + ".info.txt", 'w')
    for a in chr_tuple:
        info_file = open(sam + "." + a + ".info.temp", 'r')
        merged_info_file = open(sam + ".info.txt", 'a')
        for b in info_file:
            merged_info_file.write(b)
        info_file.close()
        os.remove(sam + "." + a + ".info.temp")
    merged_info_file.close()


if __name__ == "__main__":
    run()
    global log_file
    log_file = open(sam + "_samtopos.log", 'w')
    log_file.write(f'{time.asctime(time.localtime(time.time()))}\tlog file: Created\n')
    log_file.close()

    def append_log(module):
        log_file = open(sam + "_samtopos.log", 'a')
        log_file.write(f'{time.asctime(time.localtime(time.time()))}\t{module}: Completed\n')
        log_file.close()


    find_center(), append_log("find_center")

    def multi_processes(func):
        for i in range(0, len(chr_tuple), max_processors):
            if int(i + max_processors) > int(len(chr_tuple)):
                chr_tuple_chunk = chr_tuple[i:int(len(chr_tuple))]
            else:
                chr_tuple_chunk = chr_tuple[i:i + max_processors]
            processes = []
            for chr in chr_tuple_chunk:
                processes.append(multiprocessing.Process(target=func, args=(chr,sam, chr_name, chr_len, sigma, exclusion)))
            for process in processes:
                process.start()
            for process in processes:
                process.join()


    multi_processes(split_chr), append_log("split_chr"), os.remove(sam + ".temp")
    multi_processes(count_reads), append_log("count_reads")
    multi_processes(create_wig), append_log("create_wig")
    if rawwig:
        merge_raw_wig(), append_log("merge_raw_wig")
    multi_processes(gaussian_smoothing), append_log("gaussian_smoothing")
    multi_processes(find_peaks), append_log("find_peaks")
    if smoothwig:
        merge_smooth_wig(), append_log("merge_smooth_wig")
    multi_processes(perform_exclusion), append_log("perform_exclusion")
    multi_processes(get_info), append_log("get_info")
    merge_info(), append_log("merge_info")


    print("Done!")
