#python3 basic_aligner.py -g practice_W_1/ref_practice_W_1_chr_1.txt -r practice_W_1/reads_practice_W_1_chr_1.txt -o test_output.txt -t practice_W_1_chr_1
#python3 basic_aligner.py -g hw1_W_2/ref_hw1_W_2_chr_1.txt -r hw1_W_2/reads_hw1_W_2_chr_1.txt -o W2_output.txt -t hw1_W_2_chr_1


import sys
import argparse
import numpy as np
import time
import zipfile


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


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""
'''
kmers = []
def makekmers(read, k):
    for i in range(len(read)-k+1):
        kmers.append(read[i:i+k])
    return
'''

kmersdict = {}
def makekmers(read, k):
    for i in range(len(read)-k+1):
        if read[i:i+k] in kmersdict.keys():
            kmersdict[read[i:i+k]] += 1
        else:
            kmersdict[read[i:i+k]] = 1
    return


table = {}
def hashTable(reference, k):
    for i in range(len(reference)-k+1):
        kmer = reference[i:i+k]
        if kmer not in table.keys():
            table[kmer] = [i]
        else:
            table[kmer].append(i)   
        
def lookup(kmer):
    if kmer in table.keys():
        return table[kmer][0]
    else:
        return -1

def comparestrings(kmer, position):
    print("comparing", reference[position:position+30], "to", kmer)
    errors = 0
    minilist = []
    for i in range(30):
        if not reference[position+i] == kmer[i]:
            print("for i =", i, "and position =", position+i)
            print("reference[position+i]:", reference[position+i])
            print("kmer[i]:", kmer[i])
            errors+=1
            toadd = [reference[position+i], kmer[i], position+i]
            minilist.append(toadd)
    if (errors == 1) and not(toadd in snps):
        snps.append(toadd)

def findsnps(kmer):
    first = kmer[:10]
    second = kmer[10:20]
    third = kmer[20:]
    firstpos = lookup(first)
    secondpos = lookup(second)
    thirdpos = lookup(third)
    if (not(firstpos == -1)) and (not(thirdpos == -1)) and (firstpos+20 == thirdpos): #if first and last third is perfect match
        comparestrings(kmer, firstpos)
    if (not(firstpos == -1)) and (not(secondpos == -1)) and (firstpos+10 == secondpos):
        comparestrings(kmer, firstpos)
    if (not(secondpos == -1)) and (not(thirdpos == -1)) and (secondpos+10 == thirdpos):
        comparestrings(kmer, firstpos)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)


    """
    TODO: Call functions to do the actual read alignment here
    """
    
    #reference = reference genome
    #input_reads = list of paired reads: [['TTAGTACTGCGGTTTAAGGTAGAGCTGAGGAGGCCCCGTACAGCTTCGAA', 'TTGGGATGAAGCGGGGGTGGGCGCCGAGAGGGGCCACTATTTCCAAAACG'],...]
    snps = []
    for readpair in input_reads:
        makekmers(readpair[0], 30)
        makekmers(readpair[1], 30)

    kmers = []
    for kmer in kmersdict:
        if kmersdict[kmer] > 2:
            kmers.append(kmer)

    print(kmers)
    #print("\n\n\n\n\n\n\n")
    hashTable(reference, 10)
    print(table)
    for kmer in kmers:
        findsnps(kmer)
    #if two/three have errors, it's error
    #if one/three have error, it's SNP
    """
    TODO: Call functions to do the actual read alignment here
    """


    #snps = [['A', 'G', 3425]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
