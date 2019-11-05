from sys import exit
from sys import argv
from os import listdir as ldir
from os.path import isdir
from os.path import isfile


def main(directory):
    ncbi_identifiers = get_ncbi_identifiers(directory.rstrip('/') + '/taxonomy/ids.list')
    genomes = check_name(directory.rstrip('/') + '/genome_dir', ncbi_identifiers)
    check_fa_files(genomes, directory.rstrip('/') + '/genome_dir')
    check_blast(genomes, directory.rstrip('/') + '/blast_dir')
    check_weight(genomes, directory.rstrip('/') + '/weight_dir')


def check_name(directory, ncbi_identifiers):
    genomes = []
    ncbitogenome = {}
    naming = []
    ncbi_multi = []
    ncbi_miss = []
    for genome in ldir(directory):
        if isdir(directory + '/' + genome):
            cells = genome.split('@')
            genomes.append(genome)
            if len(cells) == 3:
                if cells[1] not in ncbitogenome:
                    ncbitogenome[cells[1]] = [genome]
                elif cells[1] in ncbi_multi:
                    ncbitogenome[cells[1]].append(genome)
                else:
                    ncbitogenome[cells[1]].append(genome)
                    ncbi_multi.append(cells[1])
            else:
                naming.append(genome)
    for ncbi_id in ncbitogenome:
        if int(ncbi_id) not in ncbi_identifiers:
            ncbi_miss.append(ncbi_id)
    print('# of genomes with wrong naming (has to be "NAME@NCBI_ID@VERSION"):\t' + str(len(naming)))
    for i in naming:
        print(i)
    print('# of NCBI ids not in taxonomy files (update taxonomy folder):\t' + str(len(ncbi_miss)))
    for i in ncbi_miss:
        print(i)
    print('# of NCBI ids that occur multiple times (remove):\t' + str(len(ncbi_multi)))
    for i in ncbi_multi:
        print(i + ':\t' + str(ncbitogenome[i]))
    return genomes


def get_ncbi_identifiers(path):
    ncbi_identifiers = set()
    with open(path, 'r') as infile:
        line = infile.readline()
        while line:
            ncbi_identifiers.add(int(line.rstrip('\n')))
            line = infile.readline()
    return ncbi_identifiers


def check_fa_files(genomes, path):
    missingfasta = []
    header = []
    for genome in genomes:
        if isfile(path + '/' + genome + '/' + genome + '.fa'):
            c1, c2 = 0, 0
            with open(path + '/' + genome + '/' + genome + '.fa', 'r') as infile:
                line = infile.readline()
                while line:
                    if line[0] == '>':
                        c1 += 1
                        if len(line) > 18:
                            c2 += 1
                    line = infile.readline()
            if c2 > 0:
                header.append((genome, c1, c2))
        else:
            missingfasta.append(genome)
    print('# of missing .fa files:\t' + str(len(missingfasta)))
    for i in missingfasta:
        print(i)
    print('# genomes with too long headers (max 18 characters):\t' + str(len(header)))
    for i in header:
        print(i[0] + ':\t' + str(i[2]) + '/' + str(i[1]))


def check_blast(genomes, path):
    missing = []
    for genome in genomes:
        if not (isfile(path + '/' + genome + '/' + genome + '.phr') and
                isfile(path + '/' + genome + '/' + genome + '.pin') and
                isfile(path + '/' + genome + '/' + genome + '.psq')):
            missing.append(genome)
    print('# of incomplete or missing blast DBs:\t' + str(len(missing)))
    for i in missing:
        print(i)


def check_weight(genomes, path):
    missing = []
    for genome in genomes:
        if not (isfile(path + '/' + genome + '/pfam.xml') and isfile(path + '/' + genome + '/cast.xml') and
                isfile(path + '/' + genome + '/tmhmm.xml') and isfile(path + '/' + genome + '/coils.xml') and
                isfile(path + '/' + genome + '/seg.xml') and isfile(path + '/' + genome + '/signalp.xml') and
                isfile(path + '/' + genome + '/smart.xml')):
            missing.append(genome)
    print('# of incomplete or missing annotations:\t' + str(len(missing)))
    for i in missing:
        print(i)


if __name__ == "__main__":
    exit(main(argv[1]))

