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
from os import listdir
from os.path import isfile, join
import time
import argparse
import subprocess
from pathlib import Path
import multiprocessing as mp
import re
from tqdm import tqdm
import h1s.h1s as h1sFn
import shutil

def main():
    version = '2.2.2'
    parser = argparse.ArgumentParser(description='You are running h1s version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', help='Input folder containing the seed sequences (protein only) in fasta format',
                            action='store', default='', required=True)
    required.add_argument('--jobName', help='Job name. This will also be file name for the output',
                            action='store', default='', required=True)
    required.add_argument('--refspec', help='Reference taxon. It should be the species the seed sequence was derived from',
                            action='store', default='', required=True)

    optional_noReusecore = parser.add_argument_group('Arguments, required when using without --reuseCore')
    optional_noReusecore.add_argument('--minDist', help='Minimum systematic distance of primer taxa for the core set compilation. Default: genus',
                            choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom'],
                            action='store', default='genus')
    optional_noReusecore.add_argument('--maxDist', help='Maximum systematic distance of primer taxa for the core set compilation. Default: kingdom',
                            choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom'],
                            action='store', default='kingdom')
    optional_noReusecore.add_argument('--coreOrth', help='Number of orthologs added to the core set. Default: 5', action='store', default=5, type=int)

    optional_paths = parser.add_argument_group('Non-default directory options')
    optional_paths.add_argument('--outpath', help='Output directory', action='store', default='')
    optional_paths.add_argument('--hmmpath', help='Path for the core ortholog directory', action='store', default='')
    optional_paths.add_argument('--blastpath', help='Path for the blastDB directory', action='store', default='')
    optional_paths.add_argument('--searchpath', help='Path for the search taxa directory', action='store', default='')
    optional_paths.add_argument('--weightpath', help='Path for the pre-calculated feature annotion directory', action='store', default='')

    addtionalIO = parser.add_argument_group('Other I/O options')
    addtionalIO.add_argument('--append', help='Append the output to existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--force', help='Overwrite existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--cleanup', help='Temporary output will be deleted. Default: True', action='store_true', default=True)
    addtionalIO.add_argument('--group', help='Allows to limit the search to a certain systematic group', action='store', default='')
    addtionalIO.add_argument('--blast', help='Determine sequence id and refspec automatically. Note, the chosen sequence id and reference species does not necessarily reflect the species the sequence was derived from.',
                                action='store_true', default=False)
    addtionalIO.add_argument('--db', help='Run oneSeq.pl in database mode. Requires a mySql database. Only for internal use.', action='store_true', default=False)

    core_options = parser.add_argument_group('Core compilation options')
    core_options.add_argument('--coreOnly', help='Compile only the core orthologs', action='store_true', default=False)
    core_options.add_argument('--reuseCore', help='Reuse existing core set of your sequence', action='store_true', default=False)
    core_options.add_argument('--coreTaxa', help='List of primer taxa that should exclusively be used for the core set compilation', action='store', default='')
    core_options.add_argument('--coreStrict', help='An ortholog is only then accepted when the reciprocity is fulfilled for each sequence in the core set',
                                action='store_true', default=False)
    core_options.add_argument('--CorecheckCoorthologsRef', help='During the core compilation, an ortholog also be accepted when its best hit in the reverse search is not the core ortholog itself, but a co-ortholog of it',
                                action='store_true', default=False)
    core_options.add_argument('--coreRep', help='Obtain only the sequence being most similar to the corresponding sequence in the core set rather than all putative co-orthologs',
                                action='store_true', default=False)
    core_options.add_argument('--coreHitLimit', help='Number of hits of the initial pHMM based search that should be evaluated via a reverse search. Default: 3',
                                action='store', default=3, type=int)
    core_options.add_argument('--distDeviation', help='The deviation in score in percent (0 = 0 percent, 1 = 100 percent) allowed for two taxa to be considered similar. Default: 0.05',
                                action='store', default=0.05, type=float)

    hamstr_options = parser.add_argument_group('Search strategy options')
    hamstr_options.add_argument('--strict', help='An ortholog is only then accepted when the reciprocity is fulfilled for each sequence in the core set',
                                action='store_true', default=False)
    hamstr_options.add_argument('--checkCoorthologsRef', help='During the final HaMStR search, accept an ortholog also when its best hit in the reverse search is not the core ortholog itself, but a co-ortholog of it',
                                action='store_true', default=False)
    hamstr_options.add_argument('--rbh', help='Requires a reciprocal best hit during the HaMStR search to accept a new ortholog',
                                action='store_true', default=False)
    hamstr_options.add_argument('--rep', help='Obtain only the sequence being most similar to the corresponding sequence in the core set rather than all putative co-orthologs',
                                action='store_true', default=False)
    hamstr_options.add_argument('--ignoreDistance', help='Ignore the distance between Taxa and to choose orthologs only based on score',
                                action='store_true', default=False)
    hamstr_options.add_argument('--lowComplexityFilterOff', help='Switch on or off the low complexity filter for the blast search. Default: False',
                                action='store_true', default=False)
    hamstr_options.add_argument('--evalBlast', help='E-value cut-off for the Blast search. Default: 0.00005',
                                action='store', default=0.00005, type=float)
    hamstr_options.add_argument('--evalHmmer', help='E-value cut-off for the HMM search. Default: 0.00005',
                                action='store', default=0.00005, type=float)
    hamstr_options.add_argument('--evalRelaxfac', help='The factor to relax the e-value cut-off (Blast search and HMM search). Default: 10',
                                action='store', default=10, type=int)
    hamstr_options.add_argument('--hitLimit', help='number of hits of the initial pHMM based search that should be evaluated via a reverse search. Default: 10',
                                action='store', default=10, type=int)
    hamstr_options.add_argument('--autoLimit', help='Invoke a lagPhase analysis on the score distribution from the hmmer search. This will determine automatically a hit limit for each query. Note, it will be effective for both the core compilation and the final ortholog search',
                                action='store_true', default=False)
    hamstr_options.add_argument('--scoreThreshold', help='Instead of setting an automatic hit limit, you can specify with this flag that only candidates with an hmm score no less than x percent of the hmm score of the best hit are further evaluated. Default: x = 10. You can change this cutoff with the option -scoreCutoff. Note, it will be effective for both the core compilation and the final ortholog search',
                                action='store', default=10, type=int)
    hamstr_options.add_argument('--scoreCutoff', help='In combination with -scoreThreshold you can define the percent range of the hmms core of the best hit up to which a candidate of the hmmsearch will be subjected for further evaluation. Default: 10',
                                action='store', default=10, type=int)
    hamstr_options.add_argument('--aligner', help='Choose between mafft-linsi or muscle for the multiple sequence alignment. DEFAULT: muscle',
                                choices=['mafft-linsi', 'muscle'], action='store', default='muscle')
    hamstr_options.add_argument('--local', help='Specify the alignment strategy during core ortholog compilation. Default: True',
                                action='store_true', default=True)
    hamstr_options.add_argument('--glocal', help='Specify the alignment strategy during core ortholog compilation. Default: False',
                                action='store_true', default=False)
    hamstr_options.add_argument('--searchTaxa', help='Specify list of search taxa', action='store', default='')

    fas_options = parser.add_argument_group('FAS options')
    fas_options.add_argument('--fasoff', help='Turn OFF FAS support', action='store_true', default=False)
    fas_options.add_argument('--countercheck', help='The FAS score will be computed in two ways', action='store_true', default=True)
    fas_options.add_argument('--coreFilter',
                                help='Specifiy mode for filtering core orthologs by FAS score. In \'relaxed\' mode candidates with insufficient FAS score will be disadvantaged. In \'strict\' mode candidates with insufficient FAS score will be deleted from the candidates list. The option \'--minScore\' specifies the cut-off of the FAS score.',
                                choices=['relaxed', 'strict'], action='store', default='')
    fas_options.add_argument('--minScore', help='Specify the threshold for coreFilter. Default: 0.75', action='store', default=0.75, type=float)

    optional = parser.add_argument_group('Other options')
    optional.add_argument('--cpu', help='Determine the number of threads to be run in parallel. Default: 4', action='store', default=4, type=int)
    optional.add_argument('--hyperthread', help='Set this flag to use hyper threading. Default: False', action='store_true', default=False)
    optional.add_argument('--showTaxa', help='Print availible taxa', action='store_true', default=False)
    optional.add_argument('--debug', help='Set this flag to obtain more detailed information about the programs actions', action='store_true', default=False)
    optional.add_argument('--silentOff', help='Show more output to terminal', action='store_true', default=False)
    optional.add_argument('--oneseqHelp', help='Print help of HaMStR-oneSeq', action='store_true', default=False)
    optional.add_argument('--oneseqVersion', help='Print version of HaMStR-oneSeq', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()

    # required arguments
    inFol = os.path.abspath(args.input)
    jobName = args.jobName
    refspec = args.refspec

    minDist = args.minDist
    maxDist = args.maxDist
    coreOrth = args.coreOrth

    # path arguments
    outpath = os.path.abspath(args.outpath)
    hmmpath = args.hmmpath
    blastpath = args.blastpath
    searchpath = args.searchpath
    weightpath = args.weightpath

    # other I/O arguments
    append = args.append
    force = args.force
    cleanup = args.cleanup
    group = args.group
    blast = args.blast
    db = args.db

    # core compilation arguments
    coreOnly = args.coreOnly
    reuseCore = args.reuseCore
    coreTaxa = args.coreTaxa
    coreStrict = args.coreStrict
    CorecheckCoorthologsRef = args.CorecheckCoorthologsRef
    coreRep = args.coreRep
    coreHitLimit = args.coreHitLimit
    distDeviation = args.distDeviation

    # hamstr arguments
    strict = args.strict
    checkCoorthologsRef = args.checkCoorthologsRef
    rbh = args.rbh
    rep = args.rep
    ignoreDistance = args.ignoreDistance
    lowComplexityFilterOff = args.lowComplexityFilterOff
    evalBlast = args.evalBlast
    evalHmmer = args.evalHmmer
    evalRelaxfac = args.evalRelaxfac
    hitLimit = args.hitLimit
    autoLimit = args.autoLimit
    scoreThreshold = args.scoreThreshold
    scoreCutoff = args.scoreCutoff
    aligner = args.aligner
    local = args.local
    glocal = args.glocal
    searchTaxa = args.searchTaxa

    # fas arguments
    fasoff = args.fasoff
    countercheck = args.countercheck
    coreFilter = args.coreFilter
    minScore = args.minScore

    # others
    cpu = args.cpu
    hyperthread = args.hyperthread
    debug = args.debug
    silentOff = args.silentOff
    if silentOff == True:
        silent = False
    else:
        silent = True
    showTaxa = args.showTaxa
    oneseqHelp = args.oneseqHelp
    oneseqVersion = args.oneseqVersion

    ### get oneSeq and data path
    oneseqPath = os.path.realpath(__file__).replace('/hms.py','')
    pathconfigFile = oneseqPath + '/bin/pathconfig.txt'
    if not os.path.exists(pathconfigFile):
        sys.exit('No pathconfig.txt found. Please run setup1s (https://github.com/BIONF/HaMStR/wiki/Installation#setup-hamstr-oneseq).')
    with open(pathconfigFile) as f:
        dataPath = f.readline().strip()
    if hmmpath == '':
        hmmpath = dataPath + '/core_orthologs'
    if blastpath == '':
        blastpath = dataPath + '/blast_dir'
    if searchpath == '':
        searchpath = dataPath + '/genome_dir'
    if weightpath == '':
        weightpath = dataPath + '/weight_dir'

    ### print oneSeq help
    if oneseqHelp:
        h1sFn.getOneseqInfo(oneseqPath, '-h')
    ### print oneSeq version
    if oneseqVersion:
        h1sFn.getOneseqInfo(oneseqPath, '-version')
    ### print available taxa
    if showTaxa:
        h1sFn.getOneseqInfo(oneseqPath, '-showTaxa')

    h1sStart = time.time()
    seeds = [f for f in listdir(inFol) if isfile(join(inFol, f))]
    print('PID ' + str(os.getpid()))
    ### run core compilation
    if reuseCore == False:
        print('Starting compiling core orthologs...')
        start = time.time()
        coreCompilationJobs = []
        for seed in seeds:
            seqFile = inFol + '/' + seed
            seqName = seed.split('.')[0]
            seqName = re.sub('[\|\.]', '_', seqName)
            ### check input arguments
            seqFile, hmmpath, blastpath, searchpath, weightpath = h1sFn.checkInput([oneseqPath, seqFile, refspec, outpath, hmmpath, blastpath, searchpath, weightpath])
            # group arguments
            basicArgs = [oneseqPath, seqFile, seqName, refspec, minDist, maxDist, coreOrth]
            ioArgs = [append, force, cleanup, group, blast, db]
            pathArgs = [outpath, hmmpath, blastpath, searchpath, weightpath]
            coreArgs = [True, reuseCore, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation]
            fasArgs = [fasoff, countercheck, coreFilter, minScore]
            hamstrArgs = [strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa]
            otherArgs = [cpu, hyperthread, debug, True]
            coreCompilationJobs.append([basicArgs, ioArgs, pathArgs, coreArgs, hamstrArgs, fasArgs, otherArgs, True])
        pool = mp.Pool(cpu)
        annoOut = []
        for _ in tqdm(pool.imap_unordered(h1sFn.h1s, coreCompilationJobs), total=len(coreCompilationJobs)):
            annoOut.append(_)
        end = time.time()
        print('==> Core compiling finished in ' + '{:5.3f}s'.format(end-start))

    ### run ortholog search
    # create list of search taxa
    print('Creating list for search taxa...')
    searchTaxa = ''
    searchGroup = 'all'
    if not group == '':
        searchTaxa = '%s/searchTaxa.txt' % (outpath)
        searchGroup = group
        cmd = 'perl %s/bin/getSearchTaxa.pl -i %s -b %s -h %s -r %s -n %s -t %s/taxonomy -o %s' % (oneseqPath, searchpath, evalBlast, evalHmmer, evalRelaxfac, searchGroup, oneseqPath, searchTaxa)
        try:
            subprocess.call([cmd], shell = True)
        except:
            sys.exit('Problem running\n%s' % (cmd))

    # do ortholog search
    mute = False
    if silent == True:
        mute = True
    if coreOnly == False:
        print('Searching orthologs for...')
        start = time.time()
        for seed in seeds:
            if mute == True:
                print(seed)
            else:
                print('\n##### ' + seed)
            seqFile = inFol + '/' + seed
            seqName = seed.split('.')[0]
            seqName = re.sub('[\|\.]', '_', seqName)
            ### check input arguments
            seqFile, hmmpath, blastpath, searchpath, weightpath = h1sFn.checkInput([oneseqPath, seqFile, refspec, outpath, hmmpath, blastpath, searchpath, weightpath])
            # group arguments
            basicArgs = [oneseqPath, seqFile, seqName, refspec, minDist, maxDist, coreOrth]
            ioArgs = [append, force, cleanup, group, blast, db]
            pathArgs = [outpath, hmmpath, blastpath, searchpath, weightpath]
            coreArgs = [False, True, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation]
            fasArgs = [False, countercheck, coreFilter, minScore]
            hamstrArgs = [strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa]
            otherArgs = [cpu, hyperthread, debug, silent]
            h1sFn.h1s([basicArgs, ioArgs, pathArgs, coreArgs, hamstrArgs, fasArgs, otherArgs, mute])
        end = time.time()
        print('==> Ortholog search finished in ' + '{:5.3f}s'.format(end-start))

    ### join output
    Path(outpath+'/'+jobName).mkdir(parents=True, exist_ok=True)
    finalFa = '%s/%s.extended.fa' % (outpath, jobName)
    with open(finalFa,'wb') as wfd:
        for seed in seeds:
            seqName = seed.split('.')[0]
            seqName = re.sub('[\|\.]', '_', seqName)
            # with open(outpath + '/' + seqName + '/' + seqName + '.extended.fa','rb') as fd:
            #     shutil.copyfileobj(fd, wfd)
            cpCmd = 'cp %s/%s/%s* %s/%s/' % (outpath, seqName, seqName, outpath, jobName)
            try:
                subprocess.call([cpCmd], shell = True)
            except:
                sys.exit('Problem running\n%s' % cpCmd)
    mergeCmd = "mergeOutput1s -i %s/%s/ -o %s" % (outpath, jobName, jobName)
    subprocess.call([mergeCmd], shell = True)
    rmCmd = 'rm -rf %s/%s' % (outpath, jobName)
    subprocess.call([rmCmd], shell = True)

    ### calculate FAS scores
    # if fasoff == False:
    #     print('Starting calculating FAS scores...')
    #     start = time.time()
    #     fasCmd = 'hamstrFAS -i %s -w %s --cores %s' % (finalFa, weightpath, cpu)
    #     try:
    #         subprocess.call([fasCmd], shell = True)
    #         end = time.time()
    #         print('==> FAS calculation finished in ' + '{:5.3f}s'.format(end-start))
    #     except:
    #         sys.exit('Problem running\n%s' % (fasCmd))


    h1sEnd = time.time()
    print('==> h1s finished in ' + '{:5.3f}s'.format(h1sEnd-h1sStart))


if __name__ == '__main__':
    main()
