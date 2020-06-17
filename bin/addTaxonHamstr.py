# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to prepare data for HaMStR oneSeq.
#  It will create a folder within genome_dir with the naming scheme of
#  HaMStR ([Species acronym]@[NCBI ID]@[Proteome version], e.g
#  HUMAN@9606@3) and a annotation file in JSON format in weight_dir
#  (optional).
#  For a long header of original FASTA sequence, only the first word
#  will be taken as the ID of new fasta file, everything after the
#  first whitespace will be removed. If this first word is not unique,
#  an automatically increasing index will be added.
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: tran@bio.uni-frankfurt.de
#
#######################################################################

import sys
import os
import argparse
from pathlib import Path
from Bio import SeqIO
import subprocess
import multiprocessing as mp
from ete3 import NCBITaxa

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def checkTaxId(taxId):
    ncbi = NCBITaxa()
    tmp = ncbi.get_rank([taxId])
    try:
        tmp = ncbi.get_rank([taxId])
        rank = tmp[int(taxId)]
        if not rank == 'species':
            print('\033[92mWARNING: rank of %s is not SPECIES (%s)\033[0m' % (taxId, rank))
    except:
        print('\033[92mWARNING: %s not found in NCBI taxonomy database!\033[0m' % taxId)


def main():
    version = '1.0.0'
    parser = argparse.ArgumentParser(description='You are running addTaxonHamstr version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-f', '--fasta', help='FASTA file of input taxon', action='store', default='', required=True)
    required.add_argument('-n', '--name', help='Acronym name of input taxon', action='store', default='', required=True, type=str)
    required.add_argument('-i', '--taxid', help='Taxonomy ID of input taxon', action='store', default='', required=True, type=int)
    required.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='', required=True)
    optional.add_argument('-v', '--verProt', help='Proteome version', action='store', default=1, type=int)
    optional.add_argument('-a', '--doAnno', help='Do annotation using annoFAS', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1', action='store', default=0, type=int)

    args = parser.parse_args()

    checkFileExist(args.fasta)
    faIn = args.fasta
    name = args.name.upper()
    taxId = str(args.taxid)
    outPath = str(Path(args.outPath).resolve())
    doAnno = args.doAnno
    ver = str(args.verProt)
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-2

    ### species name after hamstr naming scheme
    specName = name+'@'+taxId+'@'+ver

    ### create output folders
    print('Creating output folders...')
    Path(outPath + '/genome_dir').mkdir(parents = True, exist_ok = True)
    Path(outPath + '/weight_dir').mkdir(parents = True, exist_ok = True)

    ### create file in genome_dir
    print('Parsing FASTA file...')
    genomePath = outPath + '/genome_dir/' + specName
    Path(genomePath).mkdir(parents = True, exist_ok = True)
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(faIn), 'fasta')))
    f = open(genomePath + '/' + specName + '.fa', 'w')
    index = 0
    tmpDict = []
    for id in inSeq:
        if not id in tmpDict:
            tmpDict[id] = 1
        else:
            index = index + 1
            id = str(id) + '|' + str(index)
            tmpDict[id] = 1
        f.write('>%s\n%s\n' % (id, inSeq[id].seq))
    f.close()

    ### create annotation
    if doAnno:
        annoCmd = 'annoFAS -i %s/%s.fa -o %s --cpus %s' % (genomePath, specName, outPath+'/weight_dir', cpus)
        try:
            subprocess.call([annoCmd], shell = True)
        except:
            print('\033[91mProblem with running annoFAS. You can check it with this command:\n%s\033[0m' % annoCmd)

    print('Output can be found in %s within genome_dir and weight_dir folders' % outPath)

if __name__ == '__main__':
    main()
