#!/bin/env python

#######################################################################
# Copyright (C) 2019 Julian Dosch
#
# This file is part of FAS.
#
#  FAS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FAS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FAS.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


import argparse
import os

parser = argparse.ArgumentParser(description="InterPro table output parser to FAS XML input")
parser.add_argument("-i", "--input", type=str, default=None)
parser.add_argument("-o", "--output", type=str, default=None)
parser.add_argument("-s", "--singlefile", action="store_true", help="single file output")
options = parser.parse_args()



def read_input(path):
    tools = ['CDD', 'Hamap', 'PANTHER', 'Pfam', 'PIRSF', 'SUPERFAMILY', 'PRINTS', 'Gene3D', 'ProDom', 'SMART',
             'TIGRFAM', 'ProSiteProfiles', 'ProSitePatterns', 'SFLD', 'Coils', 'MobiDBLite', 'Phobius',
             'SignalP_GRAM_POSITIVE', 'SignalP_GRAM_NEGATIVE', 'SignalP_EUK', 'TMHMM']
    proteins = {}
    proteinlengths = {}
    infile = open(path, 'r')
    line = infile.readline()
    while line:
        cells = line.rstrip('\n').split('\t')
        if not cells[0] in proteins:
            proteins[cells[0]] = {}
            proteinlengths[cells[0]] = cells[2]
            for tool in tools:
                proteins[cells[0]][tool] = {}
        if not cells[4] in proteins[cells[0]][cells[3]]:
            proteins[cells[0]][cells[3]][cells[4]] = []
        proteins[cells[0]][cells[3]][cells[4]].append((cells[6], cells[7]))
        line = infile.readline()
    infile.close()
    return proteins, proteinlengths, tools


def write_output(path, proteins, proteinlengths, tools):
    if not path[-1] == '/':
        path += '/'
    if not os.path.exists(path):
        os.makedirs(path)
    for tool in tools:
        out = open(path + tool + '.xml', 'w')
        out.write('<?xml version="1.0"?>\n<tool name="' + tool + '">\n')
        for protein in proteins:
            out.write('\t<protein id="' + protein + '" length="' + proteinlengths[protein] + '">\n')
            for ftype in proteins[protein][tool]:
                out.write('\t\t<feature type="' + ftype + '" instance="' + str(len(proteins[protein][tool][ftype])) +
                          '">\n')
                for instance in proteins[protein][tool][ftype]:
                    out.write('\t\t\t<start start="' + instance[0] + '"/>\n\t\t\t<end end="' + instance[1] + '"/>\n')
                out.write('\t\t</feature>\n')
            out.write('\t</protein>\n')
        out.write('</tool>')
        out.close()


def write_output_single(path, proteins, proteinlengths):
    if not path[-1] == '/':
        path += '/'
    if not os.path.exists(path):
        os.makedirs(path)
    out = open(path + 'InterPro.xml', 'w')
    out.write('<?xml version="1.0"?>\n<tool name="InterProScan">\n')
    for protein in proteins:
        out.write('\t<protein id="' + protein + '" length="' + proteinlengths[protein] + '">\n')
        for tool in proteins[protein]:
            for ftype in proteins[protein][tool]:
                out.write('\t\t<feature type="' + ftype + '" instance="' + str(len(proteins[protein][tool][ftype])) +
                          '">\n')
                for instance in proteins[protein][tool][ftype]:
                    out.write('\t\t\t<start start="' + instance[0] + '"/>\n\t\t\t<end end="' + instance[1] + '"/>\n')
                out.write('\t\t</feature>\n')
        out.write('\t</protein>\n')
    out.write('</tool>')
    out.close()


def main(option):
    proteins, proteinlengths, tools = read_input(option.input)
    if option.singlefile:
        write_output_single(option.output, proteins, proteinlengths)
    else:
        write_output(option.output, proteins, proteinlengths, tools)

main(options)
