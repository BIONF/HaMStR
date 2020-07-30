# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to uninstall HaMStR oneSeq and its data
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
import shutil

def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)
    while True:
        # sys.stdout.write(question + prompt)
        choice = sys.stdin.readline().rstrip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with "yes" or "no" '
                             '(or "y" or "n").\n')

def main():
    version = '1.0.0'
    parser = argparse.ArgumentParser(description='You are running uninstall1s version ' + str(version) + '.')
    parser.add_argument('--data', help='Remove HaMStR-oneSeq together with all data within hamstr1s directory', action='store_true', default=False)
    args = parser.parse_args()
    data = args.data

    oneseqPath = os.path.realpath(__file__).replace('/uninstall1s.py','')
    if data:
        print('All files and folders in %s will be removed! Enter to continue' % oneseqPath)
    else:
        print('hamstr1s will be uninstalled. Its data still can be found in %s! Enter to continue' % oneseqPath)
    if query_yes_no('Are you sure?'):
        if data:
            folders = ['bin', 'core_orthologs', 'genome_dir', 'blast_dir', 'weight_dir', 'taxonomy', 'data']
            for f in folders:
                dirPath = oneseqPath+'/'+f
                if os.path.exists(os.path.abspath(dirPath)):
                    print('removing %s...' % f)
                    shutil.rmtree(dirPath)
        uninstallCmd = 'pip uninstall hamstr1s'
        subprocess.call([uninstallCmd], shell = True)
        if data:
            if os.path.exists(os.path.abspath(oneseqPath)):
                shutil.rmtree(oneseqPath)

if __name__ == '__main__':
    main()
