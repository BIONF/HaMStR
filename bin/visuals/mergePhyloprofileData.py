from sys import exit
from sys import argv
from os import listdir as ldir
from os.path import isfile


def main(directory, out):
    if not directory[-1] == '/':
        directory += '/'
    phyloprofile = open(out + '.phyloprofile', 'w')
    domains_0 = open(out + '_0.domains', 'w')
    domains_1 = open(out + '_1.domains', 'w')
    phyloprofile.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
    for infile in ldir(directory):
        if infile.endswith('.phyloprofile') and not infile == out + '.phyloprofile':
            with open(directory + infile, 'r') as reader:
                reader.readline()
                line = reader.readline()
                while line:
                    phyloprofile.write(line)
                    line = reader.readline()
        elif infile.endswith('_0.domains') and not infile == out + '_0.domains':
            with open(directory + infile, 'r') as reader:
                line = reader.readline()
                while line:
                    domains_0.write(line)
                    line = reader.readline()
        elif infile.endswith('_1.domains') and not infile == out + '_1.domains':
            with open(directory + infile, 'r') as reader:
                line = reader.readline()
                while line:
                    domains_1.write(line)
                    line = reader.readline()
    phyloprofile.close()
    domains_0.close()
    domains_1.close()


if __name__ == "__main__":
    exit(main(argv[1], argv[2]))
