from sys import exit
from sys import argv
from os import listdir as ldir
from os.path import isfile

'''
merges all Phyloprofile files (.phyloprofile, _0.domains, _1.domains) in a given directory (argument 1) as well the
extended.fa files into one file  each. The output name and directory are given in the second argument.
'''

def main(directory, out):
    if not directory[-1] == '/':
        directory += '/'
    phyloprofile = open(out + '.phyloprofile', 'w')
    domains_0 = open(out + '_forward.domains', 'w')
    domains_1 = open(out + '_reverse.domains', 'w')
    ex_fasta = open(out + '.extended.fa', 'w')
    phyloprofile.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
    for infile in ldir(directory):
        if infile.endswith('.phyloprofile') and not infile == out + '.phyloprofile':
            with open(directory + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if not line == 'geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n':
                        phyloprofile.write(line)
        elif infile.endswith('_forward.domains') and not infile == out + '_forward.domains':
            with open(directory + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    domains_0.write(line)
        elif infile.endswith('_reverse.domains') and not infile == out + '_reverse.domains':
            with open(directory + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    domains_1.write(line)
        elif infile.endswith('.extended.fa') and not infile == out + '.extended.fa':
            with open(directory + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    ex_fasta.write(line)
    phyloprofile.close()
    domains_0.close()
    domains_1.close()
    ex_fasta.close()


if __name__ == "__main__":
    exit(main(argv[1], argv[2]))
