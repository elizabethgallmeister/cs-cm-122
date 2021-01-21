#python3 basic_assembly.py -r spectrum_A_1/reads_spectrum_A_1_chr_1.txt -o spectrum_A_1_output.txt -t spectrum_A_1_chr_1
#python3 basic_assembly.py -r practice_A_2/reads_practice_A_2_chr_1.txt -o practice_A_2_output.txt -t practice_A_2_chr_1

#python3 basic_assembly.py -r reads_hw3all_A_3_chr_1.txt -o hw3all_A_3_chr_1.txt -t hw3all_A_3_chr_1

from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""

contiglist = []
def findcontig(_first):
    '''
    if out(_first) > 0:
        for oe in outpaths(_first):
            nbp = _first->oe
    '''
    #print("start kmer", _first)
    # is out(_first) > 0?
    numcycles = len(availablepaths[_first])
    #print("number of cycles is", numcycles)
    if numcycles > 0:
        for i in range(numcycles):
            #print("starting at index", i)
            contig = []
            contig.append(_first)
            #print("first is", _first)
            _next = availablepaths[_first][0]
            contig.append(_next)
            #print("next is", _next)
            availablepaths[_first].remove(_next)
            if(not availablepaths[_first]):
                del availablepaths[_first]
            while _next in onein1outnodes:
                #print("_next is a 1 in 1 out node, so we need to add stuff")
                onein1outnodes.remove(_next)
                if _next in availablepaths.keys():
                    _nn = _next
                    _next = availablepaths[_next][0]
                    availablepaths[_nn].remove(_next)
                    if(not availablepaths[_nn]):
                        del availablepaths[_nn]
                    #print("new next is", _next)
                    contig.append(_next)
                else:
                    break
                
            #print("added contig:", contig)
            contiglist.append(contig)
    #print("leaving function now")
    return

kmers = {}
def makekmers(read, k):
    for i in range(len(read)-k+1):
        if read[i:i+k] in kmers.keys():
            kmers[read[i:i+k]] += 1
        else:
            kmers[read[i:i+k]] = 1
    return

    
'''
def clean(availablepaths):
    keystodelete = []
    allkeys = availablepaths.keys()
    for kmer in allkeys:
        firsthalf = kmer[:25]
        secondhalf = kmer[25:]
        countfirst = 0
        countsecond = 0
        for key in availablepaths.keys():
            if firsthalf in key:
                countfirst += 1
            for i in range(len(availablepaths[key])):
                if firsthalf in availablepaths[key][i]:
                    countfirst += 1
            if secondhalf in key:
                countsecond += 1
            for i in range(len(availablepaths[key])):
                if secondhalf in availablepaths[key][i]:
                    countsecond += 1
        if (countfirst < 15) or (countsecond < 15):
            keystodelete.append(kmer)
    for k in keystodelete:
        del availablepaths[k]
    return availablepaths
'''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)


    availablepaths = {}
    for readpair in input_reads:
        makekmers(readpair[0], 25)
        makekmers(readpair[1], 25)
    i = 0
    for kmer in kmers:
        if kmers[kmer] > 3:
            if (kmer[:-1] in availablepaths.keys()) and (kmer[1:] not in availablepaths[kmer[:-1]]):
                availablepaths[kmer[:-1]].append(kmer[1:])
            else:
                availablepaths[kmer[:-1]] = [kmer[1:]]
        i+=1
        if i%1000 == 0:
            print(i)
    print("done")
        
    '''
    keys_sort = sorted(availablepaths.keys())
    availablepaths_sorted = {}
    for _ in keys_sort:
        availablepaths_sorted[_] = availablepaths[_]

    onein1outnodes = []
    not1in1outnodes = []
    for node in availablepaths_sorted:
        outpaths = len(availablepaths_sorted[node])
        inpaths = 0
        for others in availablepaths_sorted:
            for n in availablepaths_sorted[others]:
                if n == node:
                    inpaths = inpaths+1
        if inpaths != 1 or outpaths != 1:
            not1in1outnodes.append(node)
        else:
            onein1outnodes.append(node)
    if( len(onein1outnodes)+len(not1in1outnodes) != len(availablepaths_sorted) ):
        print("you done goofed")
    
    '''
    print(availablepaths)
    onein1outnodes = []
    not1in1outnodes = []
    for node in availablepaths:
        outpaths = len(availablepaths[node])
        inpaths = 0
        for others in availablepaths:
            for n in availablepaths[others]:
                if n == node:
                    inpaths = inpaths+1
        if inpaths != 1 or outpaths != 1:
            not1in1outnodes.append(node)
        else:
            onein1outnodes.append(node)
    if( len(onein1outnodes)+len(not1in1outnodes) != len(availablepaths) ):
        print("you done goofed")
    
    for path in not1in1outnodes:
        findcontig(path)
    for path in onein1outnodes:
        findcontig(path)

    #print(contiglist)
    contigs = []
    for minilist in contiglist:
        final = ""
        for i in range(len(minilist)-1):
            final += minilist[i][0]
        final += minilist[-1]
        contigs.append(final)

    print("\n\n\n\n\n\n\n\n\n")
    #print(contigs)
    """
            TODO: Call functions to do the actual assembly here

    """

    #contigs = ['GCTGACTAGCTAGCTACGATCGATCGATCGATCGATCGATGACTAGCTAGCTAGCGCTGACT']

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
