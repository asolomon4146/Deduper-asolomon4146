#!/usr/bin/env python

import re
import argparse
from itertools import chain

def count_lines(filename):
    with open(filename, 'r') as file:
        return sum(1 for line in file)

def get_args():
    parser = argparse.ArgumentParser(description="A program to deduplicate PCR duplicates. This PCR duplicate removal program requires a sorted SAM file of single end reads. Only the first instance of each unique read will be kept, any PCR duplicates will be omitted based on order they appear after sorting by chromosome and position.")
    parser.add_argument("-i", "--input_file", help="Your input SAM file", required=True, type = str)
    parser.add_argument("-u", "--umi_file", help="Your file contianing a list of UMIs:", required=True, type = str)
    parser.add_argument("-o", "--output_file", help="Your output file name:", required=True, type = str)
    return parser.parse_args()

args = get_args()

def desoftclipping(func_line_list: list[str], func_plus_strand: bool, plus_count: int, reverse_count: int):
    '''
    A function which takes in a SAM read split by tab into a list, and a bool stating the strandedness of the read, and returns the 5' start position of the read, accounting for soft clipping.
    '''
    #regex pattern for capturing left 5' softclipping numbers ^(\d+)S. Only 1 capture group

    func_tuple_counters_and_pos: tuple = ()
    POS_5: int = int(func_line_list[3])

    if func_plus_strand == True: # Check if plus strand
        plus_count += 1
    #Parse CIGAR string with regex for softclipping
        sc_5_match_list = re.findall(r"^(\d+)S", func_line_list[5]) #Should only match left 5'SC so only one match max.
        if sc_5_match_list:
            POS_5 = int(func_line_list[3]) - int(sc_5_match_list[0]) # True 5' start POS equals POS - left 5'SC.
    
    if func_plus_strand == False: # Check if minus strand
        reverse_count += 1
        #Parse CIGAR string with regex for softclipping
        sc_5_match_list = re.findall(r"(\d+)S$", func_line_list[5]) #Should only match right hand 5'SC so only one match max.
        M_match_list = re.findall(r"(\d+)M", func_line_list[5]) #Looking for all Ms in CIGAR
        D_match_list = re.findall(r"(\d+)D", func_line_list[5]) #Looking for all Ds in CIGAR
        #N_match_list = re.findall(r"(\d+)N", func_line_list[5]) #Looking for all Ns in CIGAR
        for x in chain(sc_5_match_list, M_match_list, D_match_list): #N_match_list):
            POS_5 +=int(x)
        #Must subtract by 1 because POS is 1 based left most mapping position not 0 based.
        POS_5 = POS_5 - 1
    
    func_tuple_counters_and_pos = (POS_5, plus_count, reverse_count)
    return func_tuple_counters_and_pos

def strand_checker(func_line_list: list[str]):
    '''
    A function which takes in a SAM read split by tab into a list, and returns a bool stating the strandedness of the read. If True, then the read is the plus strand, if False, the read is the Minus strand.
    '''
    func_plus_strand: bool = True #Initialization of strand boolean
    flag = int(func_line_list[1]) #pulls the bitwise flag out of field 2 in the SAM line and converts to int.

    if ((flag & 16) == 16) :
        func_plus_strand = False # The read is on the reverse strand (minus)
    else:
        func_plus_strand = True # The read is on the forward strand (plus)
    
    return func_plus_strand


#counting lines initially
line_count = count_lines(args.input_file)

#initialization
header_count: int = 0
Unique_read_count: int = 0
Wrong_UMI_count: int = 0
Duplicates_removed_count: int = 0
plus_counter: int = 0
minus_counter: int = 0
dupe_plus_counter: int = 0
dupe_minus_counter: int = 0

dupe_dict: dict = {}
umi_set = set({})
dupe_set = set()
current_chrom_set = set({})
current_chrom = None
c = 0
counts_dict: dict = {}

with open(args.input_file, "r") as in_fh, open(args.output_file, "w") as out_fh, open(args.umi_file, "r") as umi_fh:
    #print("Files opened")
    for i, umi_line in enumerate(umi_fh):
        umi_line = umi_line.strip()
        umi_set.add(umi_line)

    for i, line in enumerate(in_fh):
        line = line.strip()
        line_as_list = line.split()

        if line.startswith("@"):
            header_count += 1
            out_fh.write(f'{line}\n')

        else:
            current_chrom_set.add(line_as_list[2])
            plus_strand = strand_checker(line_as_list)
            POS_and_counters_tuple = desoftclipping(line_as_list, plus_strand, plus_counter, minus_counter)
            plus_counter = POS_and_counters_tuple[1]
            minus_counter = POS_and_counters_tuple[2]
            POS_5 = POS_and_counters_tuple[0]
            UMI = line_as_list[0][-8:]
            if UMI not in umi_set:
                Wrong_UMI_count += 1
            #dupe_tupe syntax: (Chromosome, UMI, strand bool, 5' SP)
            dupe_tupe: tuple = (line_as_list[2], UMI, plus_strand, POS_5)
            
            if dupe_tupe not in dupe_set and UMI in umi_set:
                if plus_strand == True:
                    dupe_plus_counter += 1
                elif plus_strand == False:
                    dupe_minus_counter += 1
                out_fh.write(f'{line}\n')

                #Populate a dictionary to hold the nuber of reads that were kept for each chromosome
                if line_as_list[2] not in counts_dict:
                    counts_dict[line_as_list[2]] = 1
                else:
                    counts_dict[line_as_list[2]] += 1

                Unique_read_count += 1
                dupe_set.add(dupe_tupe)
            elif dupe_tupe in dupe_set:
                Duplicates_removed_count += 1
            #print(f'{dupe_set}\n')

print(f'Number of header lines: {header_count}')
print(f'Number of unique reads: {Unique_read_count}')
print(f'Number of wrong UMIs: {Wrong_UMI_count}')
print(f'Number of PCR duplicates removed: {Duplicates_removed_count}')
print(f'Number of reads on the plus strand: {POS_and_counters_tuple[1]}')
print(f'Number of reads on the minus strand: {POS_and_counters_tuple[2]}')
print(f'Number of PCR duplicates on the plus strand: {int(dupe_plus_counter-POS_and_counters_tuple[1])}')
print(f'Number of PCR duplicates on the minus strand: {int(dupe_minus_counter-POS_and_counters_tuple[2])}')
for chrom in counts_dict:
    print(f'{chrom}\t{counts_dict[chrom]}')