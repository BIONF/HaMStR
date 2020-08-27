import sys
import os
from os import listdir as ldir
import argparse
from pathlib import Path

'''
merges all Phyloprofile files (.phyloprofile, _forward.domains, _reverse.domains) in a given directory (argument 1) as well the
extended.fa files into one file  each. The output name and directory are given in the second argument.
'''

def main():
    version = '1.0.0'
    parser = argparse.ArgumentParser(description='You are running mergeOutput1s version ' + str(version) + '.')
    parser.add_argument('-i','--input', help='Input directory, where all single output (.extended.fa, .phyloprofile, _forward.domains, _reverse.domains) can be found',
                        action='store', default='', required=True)
    parser.add_argument('-o','--output', help='Output name', action='store', default='', required=True)
    args = parser.parse_args()

    directory = args.input
    out = args.output
    if not os.path.exists(os.path.abspath(directory)):
        sys.exit('%s not found' % directory)
    else:
        directory = os.path.abspath(directory)

    phyloprofile = None
    domains_0 = None
    domains_1 = None
    ex_fasta = None
    for infile in ldir(directory):
        if infile.endswith('.phyloprofile') and not infile == out + '.phyloprofile':
            if not phyloprofile:
                phyloprofile = open(out + '.phyloprofile', 'w')
                phyloprofile.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if not line == 'geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n':
                        phyloprofile.write(line)
        elif infile.endswith('_forward.domains') and not infile == out + '_forward.domains':
            if not domains_0:
                domains_0 = open(out + '_forward.domains', 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    domains_0.write(line)
        elif infile.endswith('_reverse.domains') and not infile == out + '_reverse.domains':
            if not domains_1:
                domains_1 = open(out + '_reverse.domains', 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    domains_1.write(line)
        elif infile.endswith('.extended.fa') and not infile == out + '.extended.fa':
            if not ex_fasta:
                ex_fasta = open(out + '.extended.fa', 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    ex_fasta.write(line)
    if phyloprofile:
        phyloprofile.close()
    if domains_0:
        domains_0.close()
    if domains_1:
        domains_1.close()
    if ex_fasta:
        ex_fasta.close()


if __name__ == "__main__":
    sys.exit(main())
