#author: siayouyang
#for PE SAM format data
#Usage(for count reads and extract paired reads):python CountRandSelect.py -i sample.sam -p
#Usage(for random select paired reads):python CountRandSelect.py -i sample.sam -n 10000000
from optparse import OptionParser
import random
import time
import os

##OptionParser
parser = OptionParser()
parser.add_option("-i", "--input",dest="sam", help="input data in SAM format")
parser.add_option("-c", "--count", dest="count", action="store_true", help="count reads, no need to set if --paired is used")
parser.add_option("-p", "--paired", dest="paired", action="store_true", help="count reads and extract paired reads")
parser.add_option("-n", "--number", dest="number", type="int", default=0, help="number of paired reads to be selected, input must be 100% paired")
(options, args) = parser.parse_args()

sam = options.sam
count = options.count
paired = options.paired
number = options.number

def count_reads():
    time_start = time.time()
    sam_file = open(sam, 'r')
    count_file = open(sam + ".count.log", 'w')
    plus_set = set()
    minus_set = set()
    unpaired_set = set()
    paired_set = set()
    for a in sam_file:
        if a.startswith("@"):
            continue
        else:
            if int(a.split()[8]) > 0:
                plus_set.add(a.split()[0])
            elif int(a.split()[8]) < 0:
                minus_set.add(a.split()[0])
            else:
                unpaired_set.add(a.split()[0])
    unpaired_set.update(plus_set ^ minus_set)
    paired_set.update(plus_set & minus_set)
    del plus_set, minus_set
    paired_reads = len(paired_set)
    unpaired_reads = len(unpaired_set)
    count_file.write(f"statistics for {sam}:\nunpaired reads:{unpaired_reads}\npaired reads:{paired_reads}\n")
    del paired_set, unpaired_set
    sam_file.close()
    time_end = time.time()
    count_file.write(f"running time(count):{time_end - time_start}s\n")
    count_file.close()


def extract_paired():
    time_start = time.time()
    sam_file = open(sam, 'r')
    count_paired_file = open(sam + ".paired.log", 'w')
    sam_paired_file = open(sam + ".paired.sam", 'w')
    plus_set = set()
    minus_set = set()
    unpaired_set = set()
    paired_set = set()
    for a in sam_file:
        if a.startswith("@"):
            continue
        else:
            if int(a.split()[8]) > 0:
                plus_set.add(a.split()[0])
            elif int(a.split()[8]) < 0:
                minus_set.add(a.split()[0])
            else:
                unpaired_set.add(a.split()[0])
    unpaired_set.update(plus_set^minus_set)
    paired_set.update(plus_set&minus_set)
    del plus_set, minus_set
    paired_reads = len(paired_set)
    unpaired_reads = len(unpaired_set)
    count_paired_file.write(f"statistics for {sam}:\nunpaired reads:{unpaired_reads}\npaired reads:{paired_reads}\nsuccessfully extract paired reads to {sam}.paired.sam\n")
    del paired_set
    sam_file.close()
    sam_file = open(sam, 'r')
    for b in sam_file:
        if b.startswith("@"):
            sam_paired_file.write(b)
            continue
        else:
            if b.split()[0] in unpaired_set:
                continue
            else:
                sam_paired_file.write(b)
    del unpaired_set
    sam_file.close()
    sam_paired_file.close()
    time_end = time.time()
    count_paired_file.write(f"running time(extract_paired):{time_end - time_start}s\n")
    count_paired_file.close()


def random_select():
    #clean header
    sam_file = open(sam, 'r')
    sam_temp = open(sam + ".temp", 'w')
    for a in sam_file:
        if a.startswith("@"):
            continue
        else:
            sam_temp.write(a)
    sam_file.close()
    sam_temp.close()
    #sort sam
    sam_temp = open(sam + ".temp", 'r')
    sam_temp_readlines = sam_temp.readlines()
    sam_sorted = sorted(sam_temp_readlines)
    sam_sorted_temp = open(sam + ".sorted.temp", 'w')
    for b in sam_sorted:
        sam_sorted_temp.write(b)
    sam_temp.close()
    sam_sorted_temp.close()
    del sam_sorted, sam_temp_readlines
    #random select
    sam_file = open(sam, 'r')
    sam_sorted_temp = open(sam + ".sorted.temp", 'r')
    sam_sorted_temp_readlines = sam_sorted_temp.readlines()
    selected_file = open(sam + ".randsel.sam", 'w')
    for c in sam_file:
        if c.startswith("@"):
            selected_file.write(c)
            continue
        else:
            break
    line_list = []
    for d in range(0, int(len(sam_sorted_temp_readlines)/2)):
        line_list.append(int(d))
    random.shuffle(line_list)
    for e in range(0, number):
        if (number > 0) and (number <= len(line_list)):
            line = int(line_list[e])
            selected_file.write(sam_sorted_temp_readlines[int(line)*2])
            selected_file.write(sam_sorted_temp_readlines[(int(line) * 2)+1])
        else:
            print("WARNING: Inappropriate --number value(eg. exceed the number of paired reads)")
            break
    del sam_sorted_temp_readlines, line_list
    sam_file.close()
    sam_sorted_temp.close()
    selected_file.close()
    os.remove(sam + ".temp")
    os.remove(sam + ".sorted.temp")


if __name__ == "__main__":
    if count:
        count_reads()
    if paired:
        extract_paired()
    if number > 0:
        time_start = time.time()
        random_select()
        time_end = time.time()
        print(f"running time({sam}:random_select):{time_end-time_start}s\n")
    print("Done!")
