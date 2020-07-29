# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to run HaMStR oneSeq.
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
import subprocess
from pathlib import Path

def checkOptConflict(lib, conda):
    if lib:
        if (conda):
            sys.exit('*** ERROR: --lib and --conda cannot be used at the same time!')

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running runOneseq version ' + str(version) + '.')
    parser.add_argument('--conda', help='Setup HaMStR-oneSeq within a conda env', action='store_true', default=False)
    parser.add_argument('--lib', help='Install HaMStR-oneSeq libraries only', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()
    conda = args.conda
    lib = args.lib
    checkOptConflict(lib, conda)
    oneseqPath = os.path.realpath(__file__).replace('/setupOneseq.py','')
    ### run setup
    if conda:
        setupFile = oneseqPath + '/setup/setup_conda.sh'
        subprocess.call([setupFile], shell = True)
    else:
        if lib:
            setupFile = '%s/setup/setup.sh -l' % (oneseqPath)
        else:
            setupFile = '%s/setup/setup.sh' % (oneseqPath)
        subprocess.call([setupFile], shell = True)

if __name__ == '__main__':
    main()
