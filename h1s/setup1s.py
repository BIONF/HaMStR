# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to setup HaMStR oneSeq: install dependencies and
#  download pre-computed data
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
    version = '1.0.0'
    parser = argparse.ArgumentParser(description='You are running setup1s version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-o', '--outPath', help='Output path for h1s data', action='store', default='', required=True)
    optional.add_argument('--conda', help='Setup h1s within a conda env', action='store_true', default=False)
    optional.add_argument('--lib', help='Install h1s libraries only', action='store_true', default=False)
    optional.add_argument('--getPath', help='Get path to installed h1s', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()
    conda = args.conda
    lib = args.lib
    checkOptConflict(lib, conda)
    outPath = args.outPath
    Path(outPath).mkdir(parents = True, exist_ok = True)
    oneseqPath = os.path.realpath(__file__).replace('/setup1s.py','')
    ### get path
    if args.getPath:
        print(oneseqPath)
        sys.exit()
    ### run setup
    if conda:
        setupFile = '%s/setup/setup_conda.sh -o %s' % (oneseqPath, outPath)
        subprocess.call([setupFile], shell = True)
    else:
        if lib:
            setupFile = '%s/setup/setup.sh -l' % (oneseqPath)
        else:
            setupFile = '%s/setup/setup.sh -o %s' % (oneseqPath, outPath)
        subprocess.call([setupFile], shell = True)

if __name__ == '__main__':
    main()
