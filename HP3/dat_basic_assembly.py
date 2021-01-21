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
def inDegree(adj, v):
	return sum(1 if v == u else 0 for value in list(adj.values()) for u in value)
	
def outDegree(adj, v):
	return len(adj[v])

def nonBranchNode(adj, v):
	return (inDegree(adj, v) == 1) and (outDegree(adj, v) == 1)


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


    #print(input_reads) -> list of paired reads
    #[['AAAATCTGGCAATGGTCCACCATCCGAACGGTCCTGTTGCGGACCAAGGT', 'ATCAGTTAATGTTTTCCCTATTCCCTGTGTAGATAGAGGGTGGCACGGAT'],...]

    adj = defaultdict(list)
    for readpair in input_reads:
        prefix1 = readpair[0][:-1]
        prefix2 = readpair[1][:-1]
        suffix1 = readpair[0][1:]
        suffix2 = readpair[1][1:]
        adj[prefix1].append(suffix1)
        adj[prefix2].append(suffix2)
        
    contigs = []
    
    starts = [v for v in adj if not nonBranchNode(adj, v)]

    for start in starts:
        #print ('----------------------------------------------------------')
        #print ('\n\nstart a "head" node {}'.format(start))
        for v in adj[start]:
            nextV = v
            path = [start, nextV] ## take any path
            #print ('\npath to take {}'.format(path))
            #print ('is this non-branching path {}'.format(nonBranchNode(adj, nextV)))
            while nonBranchNode(adj, nextV):
                # continue to go onto the next node (node that connects to @nextV)
                # until you see a banch-path, then you exit.
                #print ('path is direct, so we continue to walk until we are able to branch')
                nextV = adj[nextV][0]
                path.append(nextV)
                #print (path)
            ##
            #print ('what is path after "while" {}'.format(path))
            r = path[0]
            r += ''.join(p[-1] for p in path[1:]) ## add string together.
            contigs.append(r)

    #
    ' '.join(sorted(contigs))
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
