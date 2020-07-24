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

def runOneseq(args):
    (f,n,i,o,c,v,a,cpus,replace,oldFAS) = args
    script = os.path.realpath(__file__).replace('addTaxaHamstr', 'addTaxonHamstr')
    cmd = 'python3 %s -f %s -n %s -i %s -o %s -v %s --cpus %s' % (script, f,n,i,o,v,cpus)
    if c == True:
        cmd = cmd + ' -c'
    if a == True:
        cmd = cmd + ' -a'
    if oldFAS == True:
        cmd = cmd + ' --oldFAS'
    if replace == True:
        cmd = cmd + ' --replace'
    # print(cmd)
    logFile = o + '/addTaxaHamstr.log'
    cmd = cmd + ' >> ' + logFile
    try:
        subprocess.call([cmd], shell = True)
    except:
        sys.exit('Problem running\n%s' % (cmd))

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running runOneseq version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', help='Input fasta file', action='store', default='', required=True)
    optional.add_argument('--setup', help='Setup HaMStR-oneSeq', action='store_true', default=False)
    optional.add_argument('--oneseqHelp', help='Print help of oneseq', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()
    input = args.input
    help = args.oneseqHelp
    setup = args.setup

    oneseqPath = os.path.realpath(__file__).replace('/runOneseq.py','')
    if help:
        if os.path.exists(oneseqPath + '/oneSeq.pl'):
            helpMsgCmd = subprocess.Popen([oneseqPath + '/hamstr.pl', '-h'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            helpMsg, helpErr = helpMsgCmd.communicate()
            print(helpMsg.decode('UTF-8'))
            print(helpErr.decode('UTF-8'))

    if setup:
        setupFile = oneseqPath + '/setup.sh'
        subprocess.call([setupFile], shell = True)


if __name__ == '__main__':
    main()
