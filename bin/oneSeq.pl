#!/usr/bin/perl

use strict;
use warnings;
use File::Copy qw(move);

use Env qw(ONESEQDIR);
use Env qw(SWATDIR);
use lib '../lib';
use Parallel::ForkManager;
use DBI;
use IO::Handle;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Term::Cap;
use POSIX;

use IPC::Run qw( run timeout );
use Time::HiRes;
use File::Path;
use File::Basename;


# Copyright (C) 2009 INGO EBERSBERGER, ebersberger@bio.uni-frankfurt.de
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 3 of the License
# or any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; If not, see http://www.gnu.org/licenses
 
# PROGRAM DESCRIPTION: oneSeq.pl is a program for targeted ortholog search in protein sequence data. 

# PROGRAM HISTORY 
## This script is based on a perl script authored by Peter Schmitzberger in the course
## of his Master's project at the CIBIV, MFPL, Vienna, Austria

## MODIFIED: 13. Aug. 2015 - solved path issues. Script will now work together with
## HaMStR 13.2.5

## Modified 14. Aug. 2015 - added the options -outpath and -hmmpath and a more refined 
## clean up after the search, in cases where a custom outpath has been chosen.

## Modified: 01. Feb. 2016: restructured major parts

## Modified: 04. Feb. 2016: - Additions - feature architecture similarity (fas) score support
##                                      - alternations in the program flow
##                                      - global/local option for alignments
##                                      - additional options
##                                      - autocleanup
##                                      - ENV SWATDIR for alignment support (local copy required or use SWATDIR=/home/holger/appz/phredphrap)
##                                      - if you run oneSeq.pl in DB mode, please adapt /bin/run-query.sh to your username and passwort
##                                      -

############ NOTES

############ General settings
my $version = 'oneSeq v.1.0.1';
##### configure
my $configure = 0;
if ($configure == 0){
	die "\n\n$version\n\nPLEASE RUN THE CONFIGURE OR CONFIGURE_MAC SCRIPT BEFORE USING oneSeq.pl\n\n";
}
##### hostname
my $hostname = `hostname`;
chomp $hostname;
#############
my $termios = new POSIX::Termios; $termios->getattr;
my $ospeed = $termios->getospeed;
my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
my ($norm, $under, $bold) = map { $t->Tputs($_,1) } qw/me md us/;
#### Paths
my $path=$ONESEQDIR;
if (!(defined $path) or !(-e $path)) {
	die "Please set the environmental variabel ONESEQDIR\n";
}  
$path =~ s/\/$//;
printDebug("Path is $path");

my $swatPath=$SWATDIR;
if (!(defined $swatPath) or !(-e $swatPath)) {
    print "You need your own copy of the alignment tool swat.\n";
    die "Please set the environmental variabel SWATDIR \n";
}
$swatPath =~ s/\/$//;
printDebug("Swat path is $swatPath"); 
#### Programs and output
my $sedprog = 'sed';
my $grepprog = 'grep';
my $alignprog = 'swat';
# my $blast_prog = 'blastall';
my $algorithm = "blastp";
my $blast_prog = 'blastp';
my $outputfmt = 'blastxml';
my $eval_blast_query = 0.0001;
my $filter = 'T';
my $annotation_prog = 'annotation.pl';
my $fas_prog = 'greedyFAS.py';
my $profile_prog = 'parseOneSeq.pl';
my $architecture_prog = 'parseArchitecture.pl';
##### ublast Baustelle: not implemented yet
my $runublast = 0;
my $ublast = 0;
my $accel = 0.8;

############ database connection details
my $dbname="123.45.678.910";	#ID adress if DB server
my $username="username";	#your user name
my $pw="password";		#your pw
my $database = "DBI:mysql:database=dbdmpng;host=$dbname";
my $getThemAll = 0;
my $updateBlast_dir = 0;

############ directory paths
my $coreOrthologsPath = "$path/core_orthologs/";
my $outputPath = "$path/output";
my $hamstrPath = "$path/bin/hamstr";
my $homeDir = $path;
my $alignmentMatrix = $path . "/BLOSUM62.txt";
my $genome_dir = "genome_dir";
my $taxaPath = "$path/$genome_dir/";
my $blastPath = "$path/blast_dir/";
my $idx_dir = "$path/taxonomy/";
my $tmpdir = "$path/tmp"; ## Baustelle
my $dataDir = $path . '/data';
my $currDir = getcwd;
my $weightPath = "$path/weight_dir/";
my $fasPath = "$path/bin/fas/";
my $visualsPath = "$path/bin/visuals/";

my @defaultRanks = ('superkingdom', 'kingdom',
        'superphylum', 'phylum', 'subphylum',
        'superclass', 'class', 'subclass', 'infraclass',
        'superorder', 'order', 'suborder', 'parvorder', 'infraorder',
        'superfamily', 'family', 'subfamily',
        'tribe', 'subtribe',
        'genus', 'subgenus',
        'species group', 'species subgroup', 'species');

#switched from online version to flatfile because it is much faster
#taxon files can be downloaded from: ftp://ftp.ncbi.nih.gov/pub/taxonomy/
my $db = Bio::DB::Taxonomy->new(-source    => 'flatfile',
                                -nodesfile => $idx_dir . 'nodes.dmp',
                                -namesfile => $idx_dir . 'names.dmp',
                                -directory => $idx_dir);
################## some variables
my $dbHandle;
my $core_hitlimit = 3; # number of hmm hits in the hamstrsearch to consider for reblast during core set generation
# number of hmm hits in the hamstrsearch to consider for reblast during final hamstr search.
# Note, this limits the number of co-orthologs that can be found.
my $hitlimit = 10;
# Setup for FAS score support (FAS support is used by default)
# Note, fas_t is set to 0.8 by default. Changes will influence sensitivity and selectivity
my $fas_support = 1;
my $countercheck = 0;
my $fasoff      = 0; 
my $fas_T       = 0.8;
my %profile     = ();
my %fas_score_keeper = ();
my $core_filter = 1;
my $eval_filter = 0.001;
my $inst_eval_filter = 0.01;

my $help;
my $showTaxa;
my $refSpec;
my $seqFile = '';
my $seqId= '';
my $seqName;
my $minDist;
my $maxDist;
my $minCoreOrthologs;
my $coreTaxa;
my $strict;
my $force;
my $group;
my $groupNode;
my $blast;
my $batch;
my $blastNode;
my $representative;
my $core_rep;
my $debug;
my $corestrict;
my $inputSeq = "";
my $rbh;
# Note, the evalue defaults ($eval_blast, $eval_hmmer) will be relaxed for final hamstr run by $eval_relaxfac
my $eval_blast = 0.00001; #1E-5
my $eval_hmmer = 0.00001; #1E-5
my $eval_relaxfac = 10; #checked in checkInput
my $coreOnly;
my $cpu = 1;
my $silent;
my $checkcoorthologsref;
my $cccr;
# Note, the alignment strategy can be local or global
# Default: local
my $local;
my $global;
my $core_filteroff;
my $dbmode = 0;         ## default run in dbmode. consider setting this in the configure step
my $vlevel = 2;         ## verbosity level
my @taxonlist = qw();
my $seqio_object;
my %taxa;
my $autoclean;
my $getversion;
my $coreex; ## flag to set when the core set already exists
my $addenv;

################# Command line options
GetOptions ("h"                 => \$help,
            "showTaxa"          => \$showTaxa,
            "refSpec=s"         => \$refSpec,
	    "db"                => \$dbmode,
	    "filter=s"          => \$filter,
            "sequence_file=s"   => \$seqFile,
            "seqId=s"           => \$seqId,
            "seqName=s"         => \$seqName,
            "silent"            => \$silent,
            "minDist=s"         => \$minDist,
            "maxDist=s"         => \$maxDist,
            "coreOrth=i"        => \$minCoreOrthologs,
            "coreTaxa=s"        => \$coreTaxa,
            "strict"            => \$strict,
            "rbh"               => \$rbh,
            "eval_blast=s"      => \$eval_blast,
            "eval_hmmer=s"      => \$eval_hmmer,
            "eval_relaxfac=s"       => \$eval_relaxfac,
            "checkCoorthologsRef"       => \$checkcoorthologsref,
            "coreCheckCoorthologsRef"   => \$cccr,
            "hitlimitHamstr=s"          => \$hitlimit,
            "coreHitlimitHamstr=s"      => \$core_hitlimit,
            "corerep"           => \$core_rep,
            "coreStrict"        => \$corestrict,
            "coreOnly"          => \$coreOnly,
            "group=s"           => \$group,
            "blast"             => \$blast,
            "batch=s"           => \$batch,
            "fas"               => \$fas_support,
            "countercheck"      => \$countercheck,
            "fasoff"            => \$fasoff,
            "corefilteroff"     => \$core_filteroff,
            "minScore=s"        => \$fas_T,
            "local"             => \$local,
            "global"            => \$global,
            "rep"               => \$representative,
            "cpu=s"             => \$cpu,
            "outpath=s"         => \$outputPath,
            "hmmpath=s"         => \$coreOrthologsPath,
            "debug"             => \$debug,
            "core_hitlimit=s"   => \$core_hitlimit,
            "hitlimit=s"        => \$hitlimit,
            "force"             => \$force,
            "cleanup"           => \$autoclean,
            "addenv=s"          => \$addenv,
            "version"           => \$getversion,
	    "reuse_core"        => \$coreex);

############# connect to the database
if ($dbmode) {
	$dbHandle = DBI->connect($database, $username, $pw)
               or die "Can not open the database!";
} 
# check additional environment
checkEnv();

############# show all taxa
if ($showTaxa) {
#get all taxa from database
#hash example: sacce_2336 -> NCBI ID for sacce_2336
   printTaxa();
   exit;
}

%taxa = getTaxa();
## debugging message
my $taxcount = keys(%taxa);
printDebug("receiving hash of taxa with $taxcount elements from sub getTaxa");
###
for (keys %taxa){
	printDebug("value of $_ is $taxa{$_}");
}
checkOptions();


my $tree = getTree();
my $treeDelFlag = 0;

if($group) {
   foreach($tree->get_nodes()) {
      if($_->id == $groupNode->id) {
         $groupNode = $_;
      }
   }
   $tree->set_root_node($groupNode);
}

if (!$coreex) {
    removeMaxDist();

    printDebug("Subroutine call removeMinDist\nRefspec is $refSpec\nTaxon is $taxa{$refSpec}\n");
    $treeDelFlag = removeMinDist($taxa{$refSpec});
}

my $outputFa =  $coreOrthologsPath . $seqName . "/" . $seqName . ".fa";
my $outputAln = $coreOrthologsPath . $seqName . "/" . $seqName . ".aln";
my $finalOutput = $dataDir . '/' . $seqName . '.extended.fa';

createFoldersAndFiles($outputFa, $seqName, $inputSeq, $refSpec);

my $curCoreOrthologs = 0;
my $hamstrSpecies = $refSpec;
my $addedTaxon = $refSpec;
my $noMoreOrthologs = 0;
my $coremode;
my %finalcontent;
my %candicontent;
# create weight_dir in oneseq's home dir (used for annotations,weighting,feature extraction)
# get annotations for seed sequence if fas support is on
if ($fas_support){
    createWeightFolder();
    getAnnotation($outputFa);
}

#core-ortholog search
if (!$coreex) {
    $coremode = 1;
    while (get_leaves($tree, $treeDelFlag) > 0 && $curCoreOrthologs < $minCoreOrthologs && $noMoreOrthologs == 0) {

        # checking the tree which determines the taxa that are going to be searched for hits
        # printDebug("Subroutine call from core-ortholog compilation\nNumber of leaves is ".get_leaves($tree)."\nCurrent core-orthologs: $curCoreOrthologs\nVar \$noMoreOrthologs is set to $noMoreOrthologs\n");
        if ($debug){
            print "\nTaxonomic Tree as text:\n";
            my $tree_as_string = $tree->as_text("tabtree");
            print $tree_as_string;
            print "\n";
        }
        
       #generate new aln
       if($curCoreOrthologs > 0) {
          createAlnMsf();
       }

       print "In round $curCoreOrthologs running hmmbuild on $outputAln\n";
       hmmbuild($coreOrthologsPath.$seqName."/hmm_dir/".$seqName.".hmm", $outputAln);

        my $pi = new Parallel::ForkManager($cpu);

        foreach my $key (get_leaves($tree)) {
            my $pid = $pi->start and next;
            print "Hamstr species: " . $key->scientific_name . " - " . @{$key->name('supplied')}[0] . "\n";
            runHamstr(@{$key->name('supplied')}[0], $seqName, $outputFa, $refSpec, $core_hitlimit, $core_rep, $corestrict, $coremode, $eval_blast, $eval_hmmer);
            $pi->finish;
        }
        $pi->wait_all_children;


       my $addedTaxon = getBestOrtholog();
       print "\n\nAdded TAXON: " . $addedTaxon . "\n\n\n\n";

       #if a new core ortholog was found
       if($addedTaxon ne "") {
            $hamstrSpecies = $hamstrSpecies . "," . $addedTaxon;

            #clear temporary result file
            if(-e $outputFa.".extended") {
                    unlink($outputFa.".extended");
            }

            #clear all alignment files
            my @files = glob("*.allscores");
            foreach my $file (@files) {
                    unlink($file);
            }
            ++$curCoreOrthologs;
            printDebug("Subroutine call from core-ortholog compilation\nTaxon is $addedTaxon\nNCBI Id is $taxa{$addedTaxon}\n");
            $treeDelFlag = removeMinDist($taxa{$addedTaxon});
       } 
       else {
            #there are no more core orthologs
            $noMoreOrthologs = 1;
            print "\nThe desired number of core orthologs could not be reached.\n";
       }
    }

    ## This is now the final round of alignment and profile hidden Markov model building
    ## It concludes the core ortholog set compilation
    if ($curCoreOrthologs < $minCoreOrthologs ){
            print "\nWARNING: The desired number of core orthologs could not be reached. Training with only $curCoreOrthologs sequences\n";
    }
    createAlnMsf();
    hmmbuild($coreOrthologsPath.$seqName."/hmm_dir/".$seqName.".hmm", $outputAln);
}
#after having calculated the core orthologous set, 
#start hamstr to find all orthologs
if (!$coreOnly) {
    $coremode = 0;
    print "Performing the final HaMStR search on all taxa\n";
    %taxa = getTaxa();
    $tree = getTree();
    #using $eval_relaxfac to relax the evalues for final hamstr
    my $final_eval_blast = $eval_blast*$eval_relaxfac;
    my $final_eval_hmmer = $eval_hmmer*$eval_relaxfac;
    if($groupNode) {
    	foreach($tree->get_nodes()) {
            if($_->id == $groupNode->id) {
                $groupNode = $_;
            }
   	}
            $tree->set_root_node($groupNode);
    }
    my $pm = new Parallel::ForkManager($cpu);

    foreach (get_leaves($tree)) {
    	my $pid = $pm->start and next;
    	runHamstr(@{$_->name('supplied')}[0], $seqName, $finalOutput, $refSpec, $hitlimit, $representative, $strict, $coremode, $final_eval_blast, $final_eval_hmmer);
    	$pm->finish;
    }
    $pm->wait_all_children;
}

## Evaluation of all orthologs that are predicted by the final hamstr run
if(!$coreOnly and $fas_support){
    print "Evaluation of predicted HaMStR orthologs.\n";
    my $processID = $$;

    ## finalOutput to hash
    open (FINORTH, "<".$finalOutput) or die "Error: Could not find $finalOutput\n";
    my $head;

    while(<FINORTH>){
        my $line = $_;
        chomp($line);
        if ($line =~ m/^>/){
            $line =~ s/>//; # clip '>' character
            $head = $line;
        }else{
            $finalcontent{$head} = $line;
        }
    }
    close (FINORTH);

    ## create array of keys
    my @k_ary;
    foreach my $k (sort keys(%finalcontent)){
        push(@k_ary, $k);
    }

    ## define chunk size for forking (Parallel::ForkManager)
    my $size = floor(scalar(@k_ary) / $cpu);

    ## create tmp folder
    my $evaluationDir = $dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/";
    mkpath($evaluationDir);

    ## handle finalcontent (final hamstr orthologs)
    ## evaluate final hamstr orthologs with FAS score
    ## $size: declares the chunk size 
    ## $evaluationDir: FAS output 
    ## @k_ary: contains all identifier (header, keys) from finalcontet (used to distribute workload to cpus)
    nFAS_score_final($size, $evaluationDir, @k_ary);

    if($autoclean){
        runAutoCleanUp($processID);
    }
    removeMetaOrthologFiles($processID);
    parseProfile($finalOutput, $seqName);
    print "XXX::".$evaluationDir."result_scores.collection".",".$finalOutput.",". $seqName.",". $dataDir."/".$seqName."\n";
    parseArchitecture($evaluationDir."result_scores.collection",$finalOutput, $seqName, $dataDir."/".$seqName);
}

## clean up the mess...
## checking for autocleanup option
if (($outputPath ne './') and !$autoclean) {
	print "\noneSeq.pl finished. Cleaning up...\n\n";
	my $answer = '';
        my $breaker = 0;
        my $del_check = 0;
	while (($answer !~ /[yn]/i) and ($breaker < 4)) {
                $breaker++;
		$answer = getInput("Do you want to remove the modified input files *.mod? [Y|N]");
                if (($breaker > 3) and ($answer !~ /[yn]/i)){
                    print "No proper answer is given: Removing modified input files *.mod\n";
                    $del_check = 1;
                }
	}
	if (($answer =~ /y/i) or ($del_check == 1)) {
		my $delcommandMod = "rm -f $outputPath/*.mod";
		system ($delcommandMod) == 0 or die "Error deleting result files\n";
	}
	$answer = '';
        $breaker = 0;
        $del_check = 0;
	while (($answer !~ /[yn]/i) and ($breaker < 4)) {
                $breaker++;
	 	$answer = getInput("Do you want to remove the tmp dir? [Y|N]");
                if (($breaker > 3) and ($answer !~ /[yn]/i)){
                    print "No proper answer is given: Removing the tmp dir\n";
                    $del_check = 1;
                }
	}
	if (($answer =~ /y/i) or ($del_check == 1)) {
		my $delcommandTmp = "rm -rf $outputPath/tmp";
		system ($delcommandTmp) == 0 or die "Error deleting result files\n";
	}
}

######################## SUBROUTINES ########################
## handle forked score calculations(sliced key array (k_ary))
## for CORE candidates
# $n: defines the size of "slices" (array with all ids is sliced into chunks)
# $e_dir: dir for FAS output file (xml)
# $c_dir: dir for candidate sequences and tmp. coreProFile
sub nFAS_score_core{

    my $n = shift;
    my $e_dir = shift;
    my $c_dir = $coreOrthologsPath . $seqName . "/fas_dir/";
    my $py = new Parallel::ForkManager($cpu);
    
    while (my @next_n = splice @_, 0, $n) {
        
        my $pid = $py->start and next;
        
        my %core_fas_0_box;

        my $ii = 0;
        while ($ii<scalar(@next_n)){
            #header: $next_n[$ii]
            #sequence: $finalcontent{$next_n[$ii]}

            my ($name,$gene_set,$gene_id,$rep_id) = split (/\|/,$next_n[$ii]); 
            my $candseqFile = $coreOrthologsPath . $seqName . "/fas_dir/" . $gene_set . "_" . $gene_id . ".candidate";

            open(CANDI_SEQ, ">".$candseqFile) or die "Error: Could not create $candseqFile\n";
            print CANDI_SEQ ">" . $next_n[$ii] . "\n" . $candicontent{$next_n[$ii]};
            close CANDI_SEQ;
            getAnnotation_Set($gene_set);

            my $headerkey = $next_n[$ii];
            # look up of annotations from the core_orthologs dir
            # running FAS with iterative Model2 (core, -o 0, seed <--vs-- hit protein)
            my $cand_annot   = getAnnotation_Candidate($gene_set,$gene_id,$candseqFile);
            my $seed_annot  = $coreOrthologsPath.$seqName."/fas_dir/annotation_dir/";
            my $mode        = 0; #FAS scoring Model2 (core)
            my $score_0 = runFAS($seed_annot.$seqName."_seed", $cand_annot, $gene_set."_".$gene_id, $seqName, $e_dir, $weightPath."/".$gene_set,$mode);
            $core_fas_0_box{$headerkey} = $score_0;
            $ii++;
        }

        ## keep child results for later
        
        keepCandidateFAS(\%core_fas_0_box, $c_dir);
        
        $py->finish;
        
    }
    $py->wait_all_children;

}
## handle forked score calculations(sliced key array (k_ary))
## for FINAL orthologs
# $n: defines the size of "slices" (array with all ids is sliced into chunks)
# $e_dir: dir for FAS output file (xml)
sub nFAS_score_final{
    my $n = shift;
    my $e_dir = shift;
    my $ps = new Parallel::ForkManager($cpu);
    
    while (my @next_n = splice @_, 0, $n) {
        
        my $pid = $ps->start and next;
        
        my %final_fas_1_box;
        my %final_fas_0_box;

        my $ii = 0;
        while ($ii<scalar(@next_n)){
            #header: $next_n[$ii]
            #sequence: $finalcontent{$next_n[$ii]}

            my ($name,$gene_set,$gene_id,$rep_id) = split (/\|/,$next_n[$ii]);            
            my $finOrth_seqFile = $e_dir . $gene_set . "_" . $gene_id . ".ortholog";

            open(ORTH_SEQ, ">".$finOrth_seqFile) or die "Error: Could not create $finOrth_seqFile\n";
            print ORTH_SEQ ">" . $next_n[$ii] . "\n" . $finalcontent{$next_n[$ii]};
            close ORTH_SEQ;
            getAnnotation_Set($gene_set);
            my $headerkey = $next_n[$ii];
            # look up of annotations from the core_orthologs dir
            # seed --vs--> hit protein
            my $cand_annot  = getAnnotation_Candidate($gene_set,$gene_id,$finOrth_seqFile);
            my $seed_annot  = $coreOrthologsPath.$seqName."/fas_dir/annotation_dir/";
            my $mode        = 1; #FAS scoring Model1 (final)
            my $score_1 = runFAS($seed_annot.$seqName."_seed", $cand_annot, $gene_set."_".$gene_id, $seqName, $e_dir, $weightPath."/".$gene_set,$mode);
            $final_fas_1_box{$headerkey} = $score_1;
            # double check the FAS results (extra calculation needed, change of direction)
            # seed <--vs-- hit protein
            if ($countercheck){
                $mode = 0;
                my $score_0 = runFAS($seed_annot.$seqName."_seed", $cand_annot, $gene_set."_".$gene_id, $seqName, $e_dir, $weightPath."/".$gene_set,$mode);
                $final_fas_0_box{$headerkey} = $score_0;
            }
            $ii++;
        }

        ## print results into profile
        printEvaluationTab(\%final_fas_1_box, \%final_fas_0_box);
        
        $ps->finish;
        
    }
    $ps->wait_all_children;
}
## parse profile
sub parseProfile{
    my ($file, $out) = ($_[0], $_[1]);
    my ($fO_base, $fO_path, $fO_suffix) = fileparse( $file, qr/\.[^.]*/ );
    my $pro_File = $fO_path."/".$fO_base.".profile";
    my $outfile = $fO_path."/".$out;
    my @cmd;
    my $pl 	= "perl";
    my $viz	= "$visualsPath/$profile_prog";
    my $i	= "-i $pro_File";
    my $o	= "-o $outfile";
    
    my ($in, $stdout, $err);
    eval {
    @cmd = ($pl,$viz,$i,$o);
    #printVariableDebug(@cmd);
    print "\n##############################\n";
    print "Writing of Visualisation file for profile.\n";

    run \@cmd, \$in, \$stdout, \$err, timeout( 600 ) or die "$profile_prog killed.\n";
    
    print "$stdout\n";
    };
    #could become debug output:
    if($err){
        print "\nERROR:\n" . $err ."\n";	
    }
}

## parse architecture from scores.collection
sub parseArchitecture{

    my ($infile, $file, $groupID, $outfile) = ($_[0], $_[1], $_[2], $_[3]);
    my ($fO_base, $fO_path, $fO_suffix) = fileparse( $file, qr/\.[^.]*/ );
    my $proFile = $fO_path."/".$fO_base.".profile";
    my @cmd;
    my $pl 	= "perl";
    my $viz	= "$visualsPath/$architecture_prog";
    my $i	= "-i=$infile";
    my $p       = "-p=$proFile";
    my $g       = "-g=$groupID";
    my $o	= "-o=$outfile";
    
    my ($in, $stdout, $err);
    eval {
    @cmd = ($pl,$viz,$i,$p,$g,$o);
    #printVariableDebug(@cmd);
    print "\n##############################\n";
    print "Writing of Visualisation file for feature architecture.\n";

    run \@cmd, \$in, \$stdout, \$err, timeout( 600 ) or die "$architecture_prog killed.\n";
    
    print "$stdout\n";
    };
    #could become debug output:
    if($err){
        print "\nERROR:\n" . $err ."\n";	
    }
}
## auto clean up can be invoked via the "-cleanup" option
# $processID: given process ID
sub runAutoCleanUp{
    my $processID = $_[0];
    print "\noneSeq.pl finished. Starting Auto Clean-up...\n\n";

    my $delCommandMod = "rm -f $outputPath/*.mod";
    system ($delCommandMod) == 0 or die "Error deleting result files\n";
    print "--> $outputPath/*.mod deleted.\n";

    my $delCommandTmp = "rm -rf $outputPath/tmp";
    system ($delCommandTmp) == 0 or die "Error deleting result files\n";
    print "--> $outputPath/tmp deleted.\n";
    
}

## removing single fasta files of predicted orthologs - meta files that are pooled in multifastas as results
# $processID: given process ID
sub removeMetaOrthologFiles{
    my $processID = $_[0];
    print "\nCompressing meta files...\n\n";

    my $delCommandCandi = "rm -f ".$coreOrthologsPath.$seqName."/fas_dir/*.candidate";
    system ($delCommandCandi) == 0 or die "Error deleting candidate files\n";
    print "--> ".$coreOrthologsPath.$seqName."/fas_dir/*.candidate\n";

    my $delCommandOrth = "rm -f ".$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/*.ortholog";
    system ($delCommandOrth) == 0 or die "Error deleting single sequence files of orthologs\n";
    print "--> ".$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/*.ortholog\n";

    ## compress *_fas.xml files in coreOrthologs
    opendir(COREFAS,$coreOrthologsPath.$seqName . "/fas_dir/fasscore_dir");
    my @fasscores = grep(/_fas\.xml/,readdir(COREFAS));
    closedir(COREFAS);
    my @content;
    my $file;
    foreach(@fasscores){
        my $catCommand = "cat ".$coreOrthologsPath.$seqName . "/fas_dir/fasscore_dir/".$_;
        push (@content, $_); 
        $file = `$catCommand`;
        push (@content, $file);
        my $delCommandXML = "rm -f ".$coreOrthologsPath.$seqName . "/fas_dir/fasscore_dir/".$_;
        system($delCommandXML) == 0 or die "Error deleting single fas score files\n";        
    }
    open(COREFASCOL,">".$coreOrthologsPath.$seqName . "/fas_dir/fasscore_dir/core_scores.collection") or die "Error: Could not create core_scores.collection\n";
    for (my $ii = 0; $ii < scalar(@content); $ii++){
        print COREFASCOL $content[$ii];
        print COREFASCOL "\n";
    }
    close(COREFASCOL);
    print "--> ".$coreOrthologsPath.$seqName."/fas_dir/fasscore_dir/*_fas.xml\n";

    ## compress *_fas.xml files in results
    opendir(COREFAS,$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir");
    @fasscores = grep(/_fas\.xml/,readdir(COREFAS));
    closedir(COREFAS);
    @content = qw();
    $file = "";
    foreach(@fasscores){
        my $catCommand = "cat ".$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/".$_;
        push (@content, $_); 
        $file = `$catCommand`;
        push (@content, $file);
        my $delCommandXML = "rm -f ".$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/".$_;
        system($delCommandXML) == 0 or die "Error deleting single fas score files\n";        
    }
    open(COREFASCOL,">".$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/result_scores.collection") or die "Error: Could not create result_scores.collection\n";
    for (my $ii = 0; $ii < scalar(@content); $ii++){
        print COREFASCOL $content[$ii];
        print COREFASCOL "\n";
    }
    close(COREFASCOL);
    print "--> ".$dataDir."/".$seqName."_".$processID."/fas_dir/fasscore_dir/*_fas.xml\n";
}
## keep FAS scores for core candidates
# %subprofile: profile (hash) of gene ids (key) and FAS scores (value) created in forked process
# $c_dir: dir (fas_dir) for core candidates and candidates.profile (tmp file) 
sub keepCandidateFAS{
    my %subprofile = %{$_[0]};
    my $c_dir = $_[1];

    my $corePro_File = $c_dir. "candidates.profile";
    # open to append
    open(PROFILE, ">>".$corePro_File) or die "Error: Could not create $corePro_File\n";
    foreach my $key (sort keys %subprofile){
        print PROFILE $key . "\t" . $subprofile{$key}. "\n";        
    }
    close PROFILE;
}

## print tab separated table of orthologs and their fas score
## writes *extended.profile in data dir
# %profile: profile (hash) of gene ids (key) and FAS scores (value, score Model 1) created in forked process
# %counterprofile: profile (hash) of gene ids (key) and FAS scores (value, score Model 2) created in forked process
sub printEvaluationTab{
    my %profile = %{$_[0]};
    my %counterprofile = %{$_[1]};
    my ($fO_base, $fO_path, $fO_suffix) = fileparse( $finalOutput, qr/\.[^.]*/ );

    my $pro_File = $fO_path."/".$fO_base.".profile";
    # open to append
    open(PROFILE, ">>".$pro_File) or die "Error: Could not create $pro_File\n";
    foreach my $key (sort keys %profile){
        print PROFILE $key . "\t" . $profile{$key};
        if ($countercheck){
            print PROFILE "\t" . $counterprofile{$key} . "\n";
        }else{
            print PROFILE "\n";
        }
    }
    close PROFILE;
}
## starting annotation_prog for given seed sequence file
# $seedseqFile: fasta file with seed sequence
sub getAnnotation {
    my ($seedseqFile) = ($_[0]);
    
    chdir($fasPath);
    my $annotationCommand = "perl $fasPath/$annotation_prog -fasta=" . $seedseqFile . " -path=" . $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir" . " -name=" . $seqName . "_seed";
    system($annotationCommand);
}
## starting annotation_prog for candidate ortholog
## annotations are saved for reuse in fas_dir/annotation_dir
## requested annotations will be extracted on demand from weight_dir/$geneset
## returns location of annotations
# $cand_geneset: gene set for taxon (candidate ortholog)
# $gene_id: gene id of candidate ortholog ($cand_geneset)
sub getAnnotation_Candidate {
    my ($cand_geneset,$gene_id,$candseqFile) = ($_[0],$_[1],$_[2]);
    my $location = '';
    
    # check for existing annotations
    # gene annotations: 
    if (-d $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir/" . $cand_geneset . "_" . $gene_id){
        # annotations already exist for candidate gene
        $location = $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir/" . $cand_geneset . "_" . $gene_id;

        return $location;
    }elsif(-d "$weightPath/$cand_geneset"){
        # annotations for $cand_geneset already exist
        #Extracting annotations from xml files in $weightPath/$cand_geneset
        chdir($fasPath);
        my $annotationCommand = "perl $fasPath/$annotation_prog -path=" . $weightPath . $cand_geneset . " -name=" . $gene_id . " -extract=" . $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir/" . $cand_geneset . "_" . $gene_id;
        system($annotationCommand);

        $location = $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir/" . $cand_geneset . "_" . $gene_id;
        return $location;
    }else{
        # no annotations exist
        print "No annotations found.\n";
        print "--> Starting to annotate gene $gene_id from species $cand_geneset.\n";
        chdir($fasPath);
        my $annotationCommand = "perl $fasPath/$annotation_prog -fasta=" . $candseqFile . " -path=" . $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir/" . " -name=" . $cand_geneset . "_" . $gene_id;
        system($annotationCommand);
        
        $location = $coreOrthologsPath . $seqName . "/fas_dir" . "/annotation_dir/" . $cand_geneset . "_" . $gene_id;
        return $location;      
    }
    
 return $location
}

## starting annotation_prog for whole gene set (weighting)
## sub getAnnotation_Set will be called if a candidate ortholog is identified
sub getAnnotation_Set{
    my ($geneset) = ($_[0]);
    chdir($fasPath);
    
    # check for existing annotations in weights_dir
    # if annotations already exist the script will skip them/no requery
    # print "Annotations for ".$geneset." will be made. This may take a while ...\n";
    my $annotationCommand = "perl $fasPath/$annotation_prog -fasta=" . $taxaPath . $geneset ."/". $geneset . ".fa -path=" . $weightPath . " -name=" . $geneset;
    system($annotationCommand);
    
}

## running actual FAS calculations via IPC
# $single: annotaions (*xml) for single protein
# $ortholog: annotations (*xml) for set (ortholog) protein
# $name: jobname for FAS
# $group: sequence name given for oneseq
# $outdir: dir for output files
# $weight: gene set of ortholog, used for weighting
# $mode: invocation mode (single --vs--> set or set --vs--> single)
sub runFAS{
    my ($single, $ortholog, $name, $group, $outdir, $weight, $mode) = ($_[0], $_[1], $_[2], $_[3], $_[4], $_[5], $_[6]);
    chdir($fasPath);

    my @cmd;
    my $py 	= "python";
    my $fas	= "$fasPath/$fas_prog";
    my $s	= "--singlepath=$single/";
    my $p	= "--proteomepath=$ortholog/";
    my $r	= "--ref_proteome=" . $weight;
    my $j	= "--jobname=$outdir/" . $name . "_". $mode ."_fas.xml";
    my $o	= "--one_vs_all=".$mode;
    my $f       = "--efilter=".$eval_filter;       #dest="efilter", default="0.001", help="E-value filter for hmm based search methods (feature based/complete sequence).")
    my $i       = "--inst_efilter=".$inst_eval_filter;  #dest="inst_efilter", default="0.01", help="E-value filter for hmm bas
    my $a	= "--raw_output=2";
    my $h	= "--help";

    my ($in, $score, $err);
    $score = "NAN";
    eval {
    @cmd = ($py,$fas,$s,$p,$r,$j,$o,$a,$f,$i);
    #printVariableDebug(@cmd);
    print "\n##############################\n";
    print "Begin of FAS score calculation.\n";
    print "--> Running ". $fas_prog ."\n";
    run \@cmd, \$in, \$score, \$err, timeout( 10000 ) or die "$fas_prog killed.\n";
    chomp($score);
    print "SCORE: $score\n";
    };
    #could become debug output:
    if($err){
        print "\nERROR:\n" . $err ."\n";	
    }
    return $score;
}

## determines the reference species and/or the sequence id of the input sequence.
sub determineRef {
    my ($infile, @refspec) = @_;
    #run blast for all available proteomes if the given sequence is not in the database   
    print "One moment please!\nLooking for the most similar sequence to your input sequence.\n\n";
    my $bestHit->{score} = 1;
    $bestHit->{evalue} = 10;
    my $outname = $$;
    ## Baustelle: Currently, we need to loop through all possible taxa to id the best matching one
    for (my $i = 0; $i < scalar(@refspec); $i++) {
        my $curTaxon = $refspec[$i];	
        ## run the blast search
        printDebug("running blast on $curTaxon");
        my $resultFile = runBlast($seqFile, $dataDir, $outname, $tmpdir, "$blastPath/$curTaxon/$curTaxon");
        my $hits = &getBestBlasthit($tmpdir, $resultFile);
        if (defined $hits and @$hits > 0) {
            #only use the best hit with the index [0]. Note, $hits is an array ref of hashrefs. 
            if($hits->[0]->{score} > $bestHit->{score}) {
                $bestHit->{name} = $hits->[0]->{name};
                $bestHit->{score} = $hits->[0]->{score};
                $bestHit->{evalue} = $hits->[0]->{evalue};
                $bestHit->{species} = $curTaxon;
            } 
        }   
    }
    return($bestHit);
}

sub checkGroup {
    my $group = shift;
    my $node = $db->get_taxon(-name => $group);
    if($node) {
        $groupNode = $node;
    } 
    else {
        print "Your selected group " . $group . " was not found in the taxonomic tree... TERMINATING\n";
        exit;
    }
}
#################################
## check in oneSeq environment for presence of required directories
## blast_dir, genome_dir, weight_dir
## changing global variables !
sub checkEnv {
    if ($addenv){
        print "Checking additional environment ".$addenv."\nPlease note that this option is in developmental status.\n";
        if ((-d "$path/weight_dir_".$addenv) and (-d "$path/genome_dir_".$addenv) and (-d "$path/blast_dir_".$addenv)){
            print "\nThe given additional environment ".$addenv." exists.\n";
        }else{
            print "\nThe given additional environment ". $addenv . " was not found...exiting.\n";
            exit;
        }
        $genome_dir = "genome_dir"."_".$addenv;
        $taxaPath = "$path/$genome_dir/";
        $blastPath = "$path/blast_dir"."_".$addenv."/";
        $weightPath = "$path/weight_dir"."_".$addenv."/";
        print "Environment set to ".$addenv."\nUsing taxa in \n\t$taxaPath\n\t$blastPath\n\t$weightPath\n\n";
    }
}
#################################
sub checkOptions {
    #### check for help
    if($help) {
        my $helpmessage = helpMessage();
        print $helpmessage;
        exit;
    }
    if ($getversion){
        print "You are running $version\n";
        print "This version supports FAS comparison.\n";
        exit;
    }
    if($eval_relaxfac < 1){
        # rethink
        if($eval_relaxfac <= 0){
            printOut("\nThe specified factor for evalue relaxation is <= 0. Please see the help text for option -eval_relaxfac. We recommend a factor > 1. Default is 10.\n",1);
            my $answer = '';
            my $breaker = 0;
            while (($answer !~ /[0-9]/i) and ($breaker < 4)) {
                $breaker++;
                $answer = getInput("Please choose a new factor (Integer) for evalue relaxation. [1,100]");
                if (($breaker > 3) and ($answer !~ /[0-9]/i)){
                    print "No proper factor is given ... exiting.\n";
                    exit;
                }
            }
            if ($answer =~ /[0-9]/i) {
                $eval_relaxfac = $answer;
            }
        }
    }
    ### check the input file
    while ((length $seqFile == 0) or ((! -e "$currDir/$seqFile") and (! -e "$dataDir/$seqFile"))) {
        if (length $seqFile > 0){
                printOut("\nThe specified file $seqFile does neither exist in current dir: $currDir or in the specified working dir $dataDir\n",1);
        }
        $seqFile = getInput("Please specify a valid file name for the seed sequence", 1);
        if ($seqFile =~ /\//){
            ## user has provided a path
            $seqFile =~ /(.*)\/(.+)/;
            if ($1) {
                ## the user has provided a relativ path
                my $relpath = $1;
                if ($relpath =~ /^\//){
                    #user has provided an absolute path
                    $dataDir = $relpath;
                }
                elsif($relpath =~ /^\.\//){
                    $dataDir = $currDir;
                }
                elsif($relpath =~ /^\.\.\//){
                    my $dataDirTmp = $currDir;
                    while ($relpath =~ /^\.\./){
                        $relpath =~ s/^\.\.\///;
                        $dataDirTmp =~ s/(.*)\/.+$/$1/;
                    }
                    $dataDir = $dataDirTmp . '/' . $relpath;
                }		
                printDebug("setting dataDir to $dataDir");
            }
            $seqFile = $2;
            printDebug("Setting infile to $seqFile in sub checkOptions");		
        }
    }
    if (-e "$currDir/$seqFile"){
        $dataDir = $currDir;
        printDebug("Setting datadir to $currDir in sub checkOptions");
    } 

    ### checking the number of core orthologs 
    while(!$minCoreOrthologs) {
        $minCoreOrthologs = getInput("Please specify the desired number of core orthologs!", 1);
        $minCoreOrthologs = checkInt($minCoreOrthologs);
    } 
    while ((!$refSpec or !$taxa{$refSpec})  && !$blast) {
        my $output = '';
        for (my $i = 0; $i < @taxonlist; $i++) {
                $output = $output . "[$i]" . "\t" . $taxonlist[$i] . "\n"; 
        }
        for (keys %taxa){
                print "value of $_ is \'$taxa{$_}\'";
        }
        printDebug("taxa contains $taxa{$refSpec}");
        my $refSpecIdx = getInput("\n" . $output . "\n" . "You did not provide a valid reference species ($refSpec). Please choose the number for the reference species your input sequence was derived from", 1);
        $refSpec = $taxonlist[$refSpecIdx];
        checkBlastDb($refSpec, $refSpec); 
    }
    ## check for blast filter
    if ($blast_prog ne 'blastall'){
        $filter = 'yes' if $filter eq 'T';
        $filter = 'no' if $filter eq 'F';
    }

    $inputSeq = fetchSequence($seqFile, $dataDir);

    ## the user has not provided a sequence id, however, the refspec is determined.
    if($seqId eq '') {
        my $besthit;
        if (!$blast){
            ## a refspec has been determined
            #use blast to search in the proteome of the specified reference species for the input sequence
            #in order to obtain a valid sequence id
            $besthit = determineRef($seqFile, ($refSpec));
        }
        else {
            $besthit = determineRef($seqFile, @taxonlist);
        }		
        $seqId = $besthit->{name};
        $refSpec = $besthit->{species};
        my $details = "Evalue: " . $besthit->{evalue};
        printOut("Seq id has been determined as $seqId in $refSpec with $details", 2);
        if($seqId eq '') {
            print "There was no significant hit for your sequence in " . $refSpec . ".\nPlease specify a sequence id on your own.\n";
            exit;
        }
    }

    if($coreTaxa) {
        if(! -e $coreTaxa) {
            print "Please specify a valid file with taxa for the core orthologs search\n";
            exit;
        }
        my @userTaxa = parseTaxaFile();
        my %newTaxa = ();
        foreach (@userTaxa) {
            $newTaxa{$_} = $taxa{$_};
        }
        $newTaxa{$refSpec} = $taxa{$refSpec};
        %taxa = %newTaxa;
    }

    if($group) {
        checkGroup($group);
    }

    if(!$seqName) {
        my $i = 0;
        while($i < 7) {
            my $j = chr(int(rand(127)));
            if($j =~ /[a-zA-Z]/) {
                $seqName .=$j;
                $i++;
            }
        }
        print "Your sequence was named: " . $seqName . "\n\n";
    }
    $outputPath = $outputPath . "/$seqName";
    if (! -d "$outputPath"){
        mkdir "$outputPath", 0777  or die "could not create the output directory $outputPath";
    }
    ## check whether a result file already exists:
    ## -force flag not set
    if (-e "$dataDir/$seqName.extended.fa" && !$force){
        my $input = '';
        my $breaker = 0;

        while (($input !~ /^[or]/i) and ($breaker < 4)) {
            $breaker++;
            $input = getInput("\nAn outputfile $dataDir/$seqName.extended.fa already exists. Shall I overwrite it [o], or rename it [r]?", 1);
            if (($breaker > 3) and ($input !~ /[or]/i)){
                    print "Please consider option -force.\n";
                    die "No proper answer is given: Quitting\n";
            }
        }
        if ($input =~ /o/i){
            unlink "$dataDir/$seqName.extended.fa" or warn "could not remove existing output file $dataDir/$seqName\n";
            printOut("Removing existing output file $dataDir/$seqName.extended.fa", 1);
        }
        else {
            move "$dataDir/$seqName.extended.fa", "$dataDir/$seqName.extended.fa.old" or warn "Could not rename existing output file $dataDir/$seqName\n";
            printOut("Renaming existing output file to $dataDir/$seqName.extended.fa.old", 2);
        }
    ## -force flag is set: existing file will be removed:
    }elsif(-e "$dataDir/$seqName.extended.fa" && $force){
        printOut ("--> $dataDir/$seqName.extended.fa already exists but will be removed due to option -force.",1);
        unlink "$dataDir/$seqName.extended.fa" or warn "Could not remove existing output file $dataDir/$seqName.extended.fa\n";
    }
    if (-e "$dataDir/$seqName.extended.profile" && !$force){
        my $input = '';
        my $breaker = 0;

        while (($input !~ /^[or]/i) and ($breaker < 4)) {
            $breaker++;
            $input = getInput("\nAn outputfile $dataDir/$seqName.extended.profile already exists. Shall I overwrite it [o], or rename it [r]?", 1);
            if (($breaker > 3) and ($input !~ /[or]/i)){
                    print "Please consider option -force.\n";
                    die "No proper answer is given: Quitting\n";
            }
        }
        if ($input =~ /o/i){
            unlink "$dataDir/$seqName.extended.profile" or warn "could not remove existing output file $dataDir/$seqName\n";
            printOut("Removing existing output file $dataDir/$seqName.extended.profile", 1);
        }
        else {
            move "$dataDir/$seqName.extended.profile", "$dataDir/$seqName.extended.profile.old" or warn "Could not rename existing output file $dataDir/$seqName\n";
            printOut("Renaming existing output file to $dataDir/$seqName.extended.profile.old", 2);
        }
    ## -force flag is set: existing file will be removed:
    }elsif(-e "$dataDir/$seqName.extended.profile" && $force){
        printOut ("--> $dataDir/$seqName.extended.profile already exists but will be removed due to option -force.", 1);
        unlink "$dataDir/$seqName.extended.profile" or warn "Could not remove existing output file $dataDir/$seqName.extended.profile\n";
    }

    my $node;   
    $node = $db->get_taxon(-taxonid => $taxa{$refSpec});
    $node->name('supplied', $refSpec);

    #### checking for the min and max distance for the core set compilation
    while (!$maxDist or (checkRank($maxDist, $node) == 0)) {
        print "You have not defined a valid maximum distance rank!\n";
        printTaxonomy($node);
        my $in = getInput('Please choose a rank by giving the number in square brackets', 1);
        $maxDist = parseInput($node, $in);
        print "You selected ". $maxDist . " as maximum rank\n\n";
    }
    while (!$minDist or (checkRank($minDist, $node) == 0)) {
        print "You have not defined a minimum distant rank!\n";
        printTaxonomy($node);
        my $in = getInput('Please choose a rank by giving the number in square brackets', 1);
        $minDist = parseInput($node, $in);
        print "You selected " . $minDist . " as minimum rank\n\n";
    }

    #### checking in fas options
    if($fasoff){
        print "You have turned FAS support off. Candidate orthologs will not be evaluated by their FAS.\n";
        # turn FAS support off
        $fas_support = 0;
    }
    ## check if user defined fas_T is off limits
    if ($fas_T < 0 or $fas_T > 1){
        print "You chose an odd fas filter (-minScore), default was 0.8.\n";
        my $answer = '';
        while ($answer < 0 or $answer > 1) {
            $answer = getInput("Please choose a fas filter [0,1] between 0 (relaxed) and 1 (stringent):");
        }
        if ($answer > 0 and $answer < 1) {
            $fas_T = $answer;
        }
    }
    ### rather strict fas filter for core orthologs: OFF
    if($core_filteroff){
        $core_filter = 0;
    }
    ### check alignment strategy
    if ($local && $global){
            print "Please specify only one alignment strategy!\n";
            print "... exiting.\n";
            exit;
    }elsif(!$local && !$global){
            print "No specific alignment strategy set. Continuing with local alignments (Smith-Waterman-Algorithm).\n";
            $local = 1;
    }

}

####################### sub check the systematic rank
sub checkRank {
	my $rank = $_[0];
	my $node = $_[1]; 
	my $rankExists = 0;   
	while($node->ancestor && $rankExists == 0) {
		if($node->rank eq $rank) {
			$rankExists = 1;
		}
		$node = $node->ancestor;
	}
	   
	if($node->rank eq $rank) {
		$rankExists = 1;
	}
	   
	return $rankExists;
}

sub createAlnMsf {
	my $linsiCommand = "linsi " . $outputFa . " > " . $outputAln;
	system($linsiCommand) == 0 or die "Could not run linsi\n";
}

################ creating folders for fas support usage
sub createWeightFolder{
    #create weight_dir in hamstr1seq home dir
    my $weightdir = $path."/"."weight_dir";
    mkdir "$weightdir", 0777 unless -d "$weightdir";
}

################
sub createFoldersAndFiles {
    my ($outputFa, $seqName, $inputSeq, $refSpec) = (@_);
    #create core orthologs directory
    my $dir = $coreOrthologsPath . $seqName;
if (!$coreex){
    mkdir "$dir", 0777 unless -d "$dir";
    my $header = $seqName . "|" . $refSpec . "|" . $seqId;
    
    #create FA file
    open (OUTPUT,  ">$outputFa") or die "Error creating fa file $outputFa\n";
    print OUTPUT  ">" . $header . "\n";
    print OUTPUT $inputSeq;
    close OUTPUT;
    
    #create the Aln file initially only with a single species in there 
    open (OUTPUT,  ">$outputAln") or die "Error creating fa file $outputAln\n";
    print OUTPUT  ">" . $header . "\n";
    print OUTPUT $inputSeq;
    close OUTPUT;	   
    
    #create the folder for the hmm output
    my $hmmdir = $dir . "/hmm_dir";
    mkdir "$hmmdir", 0777 unless -d "$hmmdir";
}
    #create the fas_dir for core orthologs if fas support is ON
    if ($fas_support){
        my $fasdir = $dir. "/fas_dir";
        mkdir "$fasdir", 0777 unless -d "$fasdir";
        
        my $annodir = $fasdir."/annotation_dir";
        mkdir "$annodir", 0777 unless -d "$annodir";
        
        my $scoredir = $fasdir."/fasscore_dir";
        mkdir "$scoredir", 0777 unless -d "$scoredir";
    }
}
#################
sub fetchSequence {
	my ($file, $filepath) = @_;
	if (! defined $filepath){
		$filepath = '.';
	}
	my $seq = "";
	open (INPUT, "<$filepath/$file") or die print "Error opening seq file\n";
	while(<INPUT>) {
		my $line = $_;
		chomp($line);
		unless($line =~ /^\>.*/) {
			$seq = $seq . $line;
		}
	}
	close INPUT;
	$seq =~ s/\s*//g;
	printOut($seq, 2);
	return $seq;
}
#################
#choose the ortholog which reaches the highest score
sub getBestOrtholog {
	printDebug("Changing to $coreOrthologsPath$seqName", "Candidate file is $outputFa".'.extended');
	chdir($coreOrthologsPath . $seqName);
	my $candidatesFile = $outputFa . ".extended";
        my %fas_box;

	if(-e $candidatesFile) {
            my $swatCommand = "$swatPath/$alignprog " . $outputFa . " " . $candidatesFile . " -M " . $alignmentMatrix . " -raw -file";
            ## adding option -nw for global alignments
            if($global){
                $swatCommand = $swatCommand. " -nw";                    
            }

            ## fas support
            if ($fas_support){

                ## candidates to hash
                open (CANDI, "<".$candidatesFile) or die "Error: Could not find $candidatesFile\n";
                my $head;
                %candicontent = ();
                while(<CANDI>){
                    my $line = $_;
                    chomp($line);
                    if ($line =~ m/^>/){
                        $line =~ s/>//; # clip '>' character
                        $head = $line;
                    }else{
                        $candicontent{$head} = $line;
                    }
                }
                close (CANDI);

                ## create array of keys
                my @k_ary;
                foreach my $k (sort keys(%candicontent)){
                    push(@k_ary, $k);
                }

                ## define chunk size for forking (Parallel::ForkManager)
                my $size = floor(scalar(@k_ary) / $cpu);

                ## create tmp folder
                my $evaluationDir = $coreOrthologsPath.$seqName."/fas_dir/annotation_dir/";
                mkpath($evaluationDir);

                ## handle candicontent (hamstr core orthologs)
                ## evaluate hamstr core orthologs with FAS score
                ## $size: declares the chunk size 
                ## $evaluationDir: FAS output for core orthologs
                ## @k_ary: contains all identifier (header, keys) from finalcontet (used to distribute workload to cpus)
                nFAS_score_core($size, $evaluationDir, @k_ary);

                # read core profile
                # to be cleared after usage
                my $coreProFile = $coreOrthologsPath.$seqName."/fas_dir/candidates.profile";
                open (CANDI, "<".$coreProFile) or die "Error: Could not find $coreProFile\n";
                while(<CANDI>){
                    my $line = $_;
                    chomp($line);
                    my @pair = split(/\t/,$line);
                    $fas_box{$pair[0]}=$pair[1];
                }
                close CANDI;
                unlink($coreProFile);
            }

            chdir($coreOrthologsPath . $seqName);
            ### Baustelle. Better to get rid of swat
            system($swatCommand) == 0 or die "Error: Could not execute swat for sequence alignment\n";
            #get all pSwat ouptut files and sum up the single scores
            #add swat scores to the fas score 
            my @files = glob("*.allscores");
            my %scores;
            my $max;
            foreach my $file (@files) {
                    open (RESULT, $file) or die "Error: Could not open file with candidate taxa\n";	 
                    $max = 0;
                    while(<RESULT>) {
                        my $line = $_;
                        chomp($line);
                        $line =~ /(.+)\s+(\d+)\s+(\d+)/;
                        my $aln_key = $1;       # first match of exp. in $line
                        my $aln_score = $3;     # thrid match of exp. in $line
                        $aln_key =~ s/\s+//g;    # clean key, no trailing whitespaces
                        chomp($aln_score);      # clean score, no trailing whitespaces

                        # keep scores and keep max score for normalization
                        if(exists $scores{$aln_key}) {
                            $scores{$aln_key} = $scores{$aln_key} + $aln_score;
                            if ($scores{$aln_key} > $max){
                                $max = $scores{$aln_key};
                            }
                        }else {
                            $scores{$aln_key} = $aln_score;
                            if ($scores{$aln_key} > $max){
                                $max = $scores{$aln_key};
                            }
                        }
                    }
                    close RESULT;
            }
            ## aln score normalization, combined score (aln/max + fas)
            ## if fas core filter is on combined score will be set to 0 if the threshold of fas_T is not met
            my $bestTaxon;
            my $bestAlnScore = 0;
            my $bestCombi_AlnFas = 0;
            my $hotCandi = 0;
            ## check for fas support
            if($fas_support){
                ## choice of best taxon: normalized aln score + fas score driven
                foreach my $key(keys%scores) { 
                    ##################################################################################print "_".$key.":::".$scores{$key}."_\n";
                    my $aln_fas;
                    if ($core_filter){
                        # set scores to 0 if they do not meet fas_T
                        if($fas_box{$key} < $fas_T){
                            #############################################################################print "_".$key.":::".$scores{$key}."_due to fas: ".$fas_box{$key}." set to 0\n";
                            $fas_box{$key}= 0;
                        }else{
                            $hotCandi++;
                        }
                        $aln_fas = ($scores{$key} / $max) + $fas_box{$key};
                        #################################################################################print $aln_fas."\n";
                    }else{
                        $aln_fas = ($scores{$key} / $max) + $fas_box{$key};
                    }
                    
                    if($aln_fas > $bestCombi_AlnFas) {
                        ##################################################################################print "set new best:_".$key.":::".$aln_fas.":::old:".$bestCombi_AlnFas."\n";
                        $bestTaxon = ">" . $key;
                        $bestCombi_AlnFas = $aln_fas;
                    }
                }
                ## Note, if no score meets fas_T the best taxon will be choosen by alignment score.
                if ($hotCandi == 0){
                    print "Warning: No candidate met the threshold $fas_T\n";
                    print "Best taxon will be chosen by alignment score only.\n";
                    $bestTaxon = chooseTaxonByAln(\%scores);
                }
            }else{
                ## choice of best taxa: alignment score driven only 
                $bestTaxon = chooseTaxonByAln(\%scores);
            }
            
            $bestTaxon =~ s/\s//g;

            open (EXTFA, $outputFa.".extended");
            my $sequenceLine = 0;
            my $bestSequence = "";

            while(<EXTFA>) {
                    my $line = $_;
                    chomp($line);
                    if($sequenceLine == 1) {
                        $bestSequence = $line;
                        $sequenceLine = -1;
                    }

                    if($line eq $bestTaxon) {
                            $sequenceLine = 1;
                    }
            }
            close EXTFA;

            my @best = split '\|', $bestTaxon;
            my $header = $best[0] . '|' . $best[1] . '|' . $best[2];
            #check for duplicates
            open (COREORTHOLOGS, ">>$outputFa") or die "Error: Could not open file: " . $outputFa . "\n";
            print COREORTHOLOGS "\n" . $header . "\n" . $bestSequence;
            close COREORTHOLOGS;
            printDebug("bestTaxon is $best[2]\n");
            ### the best taxon will be added to the primer taxon list. Create a blastdb for it
            ## Baustelle: not sure about the naming of the files here.
            checkBlastDb($best[1], $best[1]);
	    return $best[1];
	} 
	else {
		return '';
	}
}
######################
sub chooseTaxonByAln{
    my %aln_scores = %{$_[0]};
    my $bestFit;
    my $bestAlnScore = 0;
    foreach my $key(keys%aln_scores) {
        if($aln_scores{$key} > $bestAlnScore) {
            $bestFit = ">" . $key;
            $bestAlnScore = $aln_scores{$key};
        }
    }
    return $bestFit;
}
######################
sub get_leaves {
	my $tree = $_[0];
        my $delFlag = 0;
        if(defined($_[1])){
            $delFlag = $_[1];
        }
        
	my $node = $tree->get_root_node;
	my @leaves;
	my @children = ($node);
	for (@children) {
		push @children, $_->each_Descendent();
	}
	for (@children) {
		push @leaves, $_ if defined($_->name('supplied'));
	}
        # if the tree is set to be deleted
        if ($delFlag){
            @leaves = qw();
            return @leaves;
        }else{
            return @leaves;
        }
}
####### get all taxa from the database (or the $genome_dir) where a genome is available 
sub getTaxa {
	if ($dbmode) {
            my ($sql) = "select l.taxon_id, l.taxon_db, l.max_source_id, t.ncbi_id from cproteome_list.list l, taxon t where t.taxon_id = l.taxon_id and t.ncbi_id != 0";
            my ($query) = $dbHandle->prepare($sql);
            $query->execute();
            while(my @result = $query->fetchrow_array) {
                ## modified by ingo: make sure to capture the max_source_id
                my $tax_src = $result[1] . '@' . $result[3] . '@' . $result[2];
                push @taxonlist, $tax_src;
                $taxa{$tax_src} = $result[3];
                printDebug("ncbiid of $tax_src is $taxa{$tax_src}");
                if ($getThemAll){
                    getProteome($tax_src);
                }
            }
	}
	else {
		@taxonlist = `ls $path/$genome_dir`;
		chomp @taxonlist;
		for (my $i = 0; $i < @taxonlist; $i++) {
			my ($taxon_name, $ncbi_id, $src_id) = split /@/, $taxonlist[$i];
			if (!$src_id) {
				$src_id = '';
			}
			$taxon_name = $taxonlist[$i];
			$taxa{$taxon_name} = $ncbi_id;
		}
	}
	### if the blast option is chosen, we will need blast databases for all taxa
	### Baustelle: have one database including all taxa to run just a single instead of n blast searches
	if ($blast or $updateBlast_dir){
		for (my $i = 0; $i < @taxonlist; $i++){
			checkBlastDb($taxonlist[$i], $taxonlist[$i]);
		}
                if ($updateBlast_dir){
                    print "\nMissing blast databases updated. Exiting.\n";
                    exit;
                }
	}
	my $hashcount = keys(%taxa);
	printDebug("Returning $hashcount taxa from subroutine getTaxa");
	return(%taxa);
}
####################
sub getTree {
	# the full lineages of the species are merged into a single tree
	my $tree;
	foreach my $key (keys%taxa) {
		my $node = $db->get_taxon(-taxonid => $taxa{$key});
                print "key: ".$key."\n";
		$node->name('supplied', $key);
		if($tree) {
			$tree->merge_lineage($node);
		} 
		else {
			$tree = Bio::Tree::Tree->new(-verbose => $db->verbose, -node => $node);
		}
	}
        if ($debug){
            print "\nTaxonomic Tree as text:\n";
            my $tree_as_string = $tree->as_text("tabtree");
            print $tree_as_string;
            print "\n";
        }
	return $tree;
}


##################### perform the hamstr search for orthologs
# using the core-orthologs found in the previous steps
sub runHamstr {
	my ($taxon, $seqName, $outputFa, $refSpec, $hitlimit, $rep, $sub_strict, $subcoremode, $ev_blst, $ev_hmm) = (@_);
	my $taxaDir = $taxaPath . $taxon;
	printDebug("Running hamstr: $taxon\t$seqName\t$outputFa\t$refSpec\t$taxaDir"); 
	if (! -e $taxaDir) {
		## backward compatibility. I used to name the dirs with the ending .dir
		if (-e "$taxaDir.dir"){
			$taxaDir = $taxaDir . '.dir';
		}
	} 
	$taxaDir =~ s/\s*//g;
	if(! -e $taxaDir and $dbmode) {
		getProteome($taxon);
	}
	if (-e $taxaDir){
		print "hamstr for taxon: " . $taxon . "\n";
		chdir($taxaDir) or die "Error: Directory for " . $taxon  . " does not exist!\n";
		my $seqfile = $taxon . ".fa";

		if(! -e $seqfile) {
			printOut("Could not find $seqfile. Check naming conventions of the files. Exiting...");
			exit;
		}
	      
		if($seqFile ne "") {
		    my $taxon_id = substr($taxon, 6, length($taxon));
		    my @hamstr = ($hamstrPath, "-sequence_file=".$seqfile, "-fasta_file=".$outputFa, "-hmmpath=".$coreOrthologsPath , "-outpath=".$outputPath,
			     "-blastpath=".$blastPath , "-protein", "-hmmset=".$seqName, "-taxon=".$taxon, "-force", 
			     "-eval_blast=".$ev_blst, "-eval_hmmer=".$ev_hmm, "-central");

		    my $resultFile;
                    if (defined $hitlimit) {
                        push(@hamstr, "-hit_limit=$hitlimit"); 
                    }
		    if($sub_strict) {
		    	push(@hamstr, "-strict");
		    	$resultFile = $outputPath . "/fa_dir_" . $taxon . '_' . $seqName . "_strict/" . $seqName . ".fa";
		    }
		    else {
			push(@hamstr,  "-refspec=".$refSpec);
			$resultFile = $outputPath . "/fa_dir_" . $taxon . '_' . $seqName . "_" . $refSpec . "/" . $seqName . ".fa";
		    }
		    if($rep) {
			push(@hamstr, "-representative");
		    }
                    if ($checkcoorthologsref and $subcoremode==0){
                        push @hamstr, '-checkCoorthologsRef';
		    }
                    if ($cccr and $subcoremode==1){
                        push @hamstr, '-checkCoorthologsRef';
		    }
		    if ($rbh) {
                        push @hamstr, "-rbh";
		    }
			
		    if ($silent) {
                        push @hamstr, "-silent";
		    }
                    printDebug(@hamstr);
	 
		    system(@hamstr) == 0 or die "Error: hamstr failed for " . $taxon . "\n";

		    if(-e $resultFile) {
		    	if ($outputFa !~ /extended/){
                            $outputFa .= '.extended';
			}
		    	## Baustelle: check that this also works with the original hamstrcore module as here a tail command was used.
			my $tailCommand = "grep -A 1 '$taxon.*|[01]\$' " . $resultFile . "|sed -e 's/\\([^|]\\{1,\\}\\)|[^|]*|\\([^|]\\{1,\\}\\)|\\([^|]\\{1,\\}\\)|\\([01]\\)\$/\\1|\\2|\\3|\\4/' >>" . $outputFa;
			printDebug("Post-processing of HaMStR\n$tailCommand\n");
			system($tailCommand);
		    }
		    else {
			printDebug("$resultFile not found");
		    }
		}
		#remove the created folders and files
		#delete fa_dir
		my $delCommandFa;
		my $delCommandHmm;
		my $delCommandHam;
		
		if (!$strict) {
			$delCommandFa = "rm -rf " . $outputPath . "/fa_dir_" . $taxon . "_" . $seqName . "_" . $refSpec;
			$delCommandHmm = "rm -rf " .  $outputPath . "/hmm_search_" . $taxon . "_" . $seqName;
			$delCommandHam = "rm -f " . $outputPath . "/hamstrsearch_" . $taxon . "_" . $seqName . ".out";
		}
		else {
			$delCommandFa = "rm -rf " . $outputPath . "/fa_dir_" . $taxon . "_" . $seqName . "_strict";
			$delCommandHmm = "rm -rf " .  $outputPath . "/hmm_search_" . $taxon . "_" . $seqName;
			$delCommandHam = "rm -f " . $outputPath . "/hamstrsearch_" . $taxon . "_" . $seqName . ".strict.out";
		}
		printDebug("executing $delCommandFa", "executing $delCommandHmm", "executing $delCommandHam");
		if (!$debug) {
			system ($delCommandFa) == 0 or die "Error deleting result files\n";
			system ($delCommandHmm) == 0 or die "Error deleting result files\n";  
			system ($delCommandHam) == 0 or die "Error deleting result files\n";
		}
	}
	else {
		print "No protein set available for $taxon. Failed to fetch it from database and nothing at $taxaDir. Skipping!\n";
	}
}

##########################
sub hmmbuild {
	my @hmmbuild = ("hmmbuild", $_[0], $_[1]);
	system(@hmmbuild) == 0 or die "hmmbuild failed";
}

sub parseInput {
	my $node =  $_[0];
	my $level = $_[1];
	my $rank = $node->rank;
	#print "\nLEVEL:".$level."\n";
        #print "\nRANK:".$rank."\n";    
	while($level > 0) {
               #print "\nLEVEL:".$level."\n";
		$node = $node->ancestor;
		#if($node->rank ne "no rank") {
			$rank = $node->rank;
                        --$level;
		#}
	}
        print "\nRETURN RANK: ".$rank."\n";
	return $rank;
}
##########################
sub parseTaxaFile {
	open (INPUT, "<$coreTaxa") or die print "Error opening file with taxa for core orthologs search\n";
	my @userTaxa;
	while(<INPUT>) {
		my $line = $_;
		chomp($line);
		if(!$taxa{$line}) {
			print "You specified " . $line . " in your core orthologs file but the taxon is not in the database!\n";
			exit;
		} 
		else {
			push(@userTaxa, $line);
		}
	}
	close INPUT;
	return @userTaxa;
}
##########################
sub printTaxa {
	my @result = qw();
	if ($dbmode) {
		print "taxon_schema\tsource_id\ttaxon name\n";
		print "------------\t---------\t----------\n";
		my ($sql) = "select t.name, c.taxon_db, c.max_source_id from taxon t, cproteome_list.list c where t.taxon_id=c.taxon_id";
		my ($query) = $dbHandle->prepare($sql);
		$query->execute();
		@result = $query->fetchrow_array;
		while(my @result = $query->fetchrow_array) {
			print $result[1] . " \t" . $result[2] . "\t" . $result[0] . "\n";
		}
	}
	else {
		print "Taxon_Name\tNCBI_ID\n";
		print "-------------\t------------\n";
                my $taxacall= "ls $path/$genome_dir |$sedprog -e 's/@/\t/'";
		@result = `$taxacall`;
		chomp @result;
		print join "\n", @result;
		print "\n";
	}
}
###########################
sub printTaxonomy {
	my $node = $_[0];
	my $i = 0;
	if($node->rank eq "species") {
		print "[" . $i . "]: " . $node->rank . " (" . $node->scientific_name . ")\n";
		while($node->ancestor) {
			$node = $node->ancestor;
			#if($node->rank ne "no rank") {
				++$i;
				print "[" . $i . "]: " . $node->rank . " (" . $node->scientific_name . ")\n";
			#}
		}      
	}
}
############################
sub remove_branch {
	my $node = $_[0];
        my $delFlag = 0;
	printDebug("Subroutine remove_branch\nNode is $node\nRank of node: ".$node->rank."\nNumber of leaves before removing branch ".get_leaves($tree)."\n\n");   
        
        # undef the tree if there is only one leave left which must be removed
        if (get_leaves($tree) == 1){
            $delFlag = 1;
        }else{
            while (defined $node->ancestor) {
                    last if $node->ancestor->each_Descendent > 1;
                    $node = $node->ancestor;
            }
            $node->remove_all_Descendents;
            if(defined $node->ancestor) {
                    $node->ancestor->remove_Descendent($node);
            }
        }
        printDebug("Subroutine remove_branch\nNode is $node\nRank of node: ".$node->rank."\nNumber of leaves after removing branch ".get_leaves($tree, $delFlag)."\n\n");   
        return $delFlag;
}
############################
sub removeMaxDist {
	my $node = $tree->find_node(-ncbi_taxid => $taxa{$refSpec});
	my $root = $tree->get_root_node();
        
        if ($maxDist eq "no rank"){
            $tree->set_root_node($root);
        }else{
            while($node->rank ne $maxDist && $node != $root) {
                    $node = $node->ancestor;
            }
            $tree->set_root_node($node);
        }
}
############################
# node determines the node in the tree in accordance to the given ncbi taxon id
sub removeMinDist {
	my $ncbiId = $_[0];
	my $node = $tree->find_node(-ncbi_taxid => $ncbiId);
	my $root = $tree->get_root_node();
        my $delFlag;

	printDebug("Subroutine removeMinDist\nncbiID is $ncbiId\nNode is  $node\nRank of node is ".$node->rank."\nroot is $root\nMinimal distance is $minDist\n");
		
        # increasing the rank of the node
	while($node->rank ne $minDist && $node != $root && defined($node->ancestor)) {
            if ($debug){
                print "Increasig the rank\nRank: ".$node->rank."\nNode: ".$node."\n\n";
            }   

            $node = $node->ancestor;             
	}
	   
	#if the species has the same ranks as the references species
	if($node == $root) {
		my @toCompare = ();
		my $i = @defaultRanks - 1;
		while($i >= 0 && $defaultRanks[$i] ne $minDist) {
			push(@toCompare, $defaultRanks[$i]);
			--$i;
		}
		$node = $tree->find_node(-ncbi_taxid => $ncbiId);
		my $lastToCompare = $toCompare[$#toCompare];
		foreach(@toCompare) {
			while($node->rank eq "no rank") {
			$node = $node->ancestor;
			}
			if($node->rank ne $lastToCompare && $node->rank eq $_) {
				$node = $node->ancestor;
			}
		}
	}
	$delFlag = remove_branch($node);
        return $delFlag;
}
##########################
sub getProteome {
	my $taxstring = shift;
	my $outdir = $taxstring;
	my $outfile = $taxstring;
	my @outfile;
	$taxstring =~ /(.*)@([0-9]+)@([0-9]+)/;
	my ($schema, $ncbi_id, $src_id) = ($1, $2, $3);
	print "\n\nAttempting to fetch information for $schema using source id $src_id\n\n";
		
	## create the relevant directory
	if (!-e "$taxaPath/$outdir"){
		print "creating directory $taxaPath/$outdir\n";
		mkpath($taxaPath."/".$outdir);
		if (-e "$taxaPath/$outdir") {
			print "succeeded\n";
		}
		else {
			print "create directory failed\n";
		}
	}
	#    
	### This is the sql statement required for fetching the sequence information from the database
        ## Using Here Documents #######
        my $sql = <<"________END_OF_STATEMENT";
        use $schema;
        select concat('>',i.id, '\\n', p.seq) from ids i, protein p
        where
        i.protein_id = p.id and
        i.representative = 1 and
        i.src_id = $src_id and
        length(p.seq) > 30;
________END_OF_STATEMENT
        ## the previous line must be exactly like it is to match the end note (Here Doc)

	printDebug("$sql\n");
	open (OUTQUERY, ">$taxaPath/$outdir/query.sql") or die "Could neither find nor create query.sql in $taxaPath/$outdir";
	print OUTQUERY $sql;
	close OUTQUERY;	    
	print "attempting to enter $taxaPath/$outdir\n";
	chdir("$taxaPath/$outdir") or die "could not enter $taxaPath/$outdir";
	`$homeDir/bin/run-query.sh $schema $ncbi_id $src_id`;
}
############ 
#Baustelle: run generation of BlastDb in a sub routine
	## now create the relevant blast directories if necessary
sub checkBlastDb {
	my ($taxstring, $filename) = @_;
	## $taxstring identifies the species directory, $filename identifies the name of the file containing the protein set
	while (! -e "$taxaPath/$taxstring/$filename.fa"){
		my $count = 0; ## avoid endless loop
		printDebug("could not find $taxaPath/$taxstring/$filename.fa\n");
		getProteome($taxstring);
		if ($count == 5){
			die "could not find $taxaPath/$taxstring/$filename.fa and could not retrieve this information from the database.\nTerminating...\n\n";
		}	
	}
	if (! -e "$blastPath/$taxstring" or $updateBlast_dir){
		`mkdir $blastPath/$taxstring`;
	}
	if (! -e "$blastPath/$taxstring/$filename.fa" or $updateBlast_dir){
		`ln -s $taxaPath/$taxstring/$filename.fa $blastPath/$taxstring/$filename.fa`;
	}
	if (! -e "$blastPath/$taxstring/$filename.pin" or $updateBlast_dir){
		chdir("$blastPath/$taxstring") or die "failed to change to dir\n";
		if ($blast_prog eq 'blastall'){
			`formatdb -i $filename.fa -t $filename -n $filename`;
		}
		elsif ($blast_prog eq 'blastp'){
			printOut("attempting to run makeblastdb", 2);
			`makeblastdb -in $filename.fa -dbtype prot -title $filename -out $filename`;
		}		
	}
}
#################
sub printDebug{
	my @message = @_;
	if ($debug){
		print join "\n", @message;
		print "\n";
	}
}
sub printVariableDebug{
    my @values = @_;
    print "\n\nDEBUG\n";
    foreach (@values){
        print $_."\n";
    }
    print "\nEND OF DEBUG\n\n";
}
#################
sub getInput {
	my ($message, $dieopt) = @_;
	if ($dieopt){
		$message .= ', or type \'q\' to quit'; 
	}
	print ("\n" . $message . ": ");
	my $input = <STDIN>;
	chomp $input;
	if ($input =~ /^q$/i and $dieopt) {
		die "Quitting!\n";
	}
	else {
		return ($input);
	}
} 
#################
sub runBlast {
	my ($query, $inpath, $outname, $outpath, $blastdb) = @_;
			printDebug("running $blast_prog on database $blastdb using input $inpath/$query and writing to $outpath/$outname.blast");

	if ($blast_prog =~ /blast[px]/) {
		!`$blast_prog -db $blastdb -seg $filter -max_target_seqs 10 -evalue $eval_blast_query -outfmt 5 -query $inpath/$query -out $outpath/$outname.blast` or die "Problem running $blast_prog\n";	
	}
	elsif ($blast_prog =~ /blastall/) {
		!`$blast_prog -p $algorithm -d $blastdb -F $filter -e $eval_blast_query -m7 -i $inpath/$query -o $outpath/$outname.blast` or die "Problem running $blast_prog\n"
	}
	else {
		`$blast_prog -ublast $inpath/$query -db $blastdb -accel $accel -evalue $eval_blast_query -blast6out $outpath/$outname.blast` or die "Problem running $blast_prog\n"; 
			
		## sort the output as ublast does not do it (at least not for ESTs)
		`sort -n -r -k 12 $outpath/$outname.blast >$outpath/blastsort.tmp`;
		`mv $outpath/blastsort.tmp $outpath/$outname.blast`;
	}
	printDebug("returning $outname.blast for subroutine runBlast\n");
	return("$outname.blast");
}
#############
sub getBestBlasthit {
	my $hits;
	my $frame;
	my $count = 0;
	my ($inpath, $resultfile) = @_;
	printDebug("Sub getBestBlasthit running on $inpath/$resultfile");
	my $searchio = Bio::SearchIO->new(
	-file        => "$inpath/$resultfile",
	-format      => $outputfmt) 
	or die "parse failed";
	while(my $result = $searchio->next_result){
		my $sig;
		my $sig_old;
		while( my $hit = $result->next_hit) {
			my $frameval = $hit->strand('query');
			if ($frameval >0){
				$frame = '+';
			}
			elsif ($frameval <0 ) {
				$frame = '-';
			}
			else {
				$frame = 'na';
			}
			## now I enter all top hits having the same score into the result
			$sig = $hit->score;	
			if (!defined $sig_old) {
				$sig_old = $sig;
			}
			if ($sig == $sig_old) {
				$hits->[$count]->{name} = $hit->name;
				$hits->[$count]->{score} = $sig;
				$hits->[$count]->{evalue} = $hit->significance;
				$count ++;
			}
			else {
			## there is no lower ranking hit with the same score as the best hit. End the loop.
			last;
	  	  }
		}
    }
    return($hits);
}
##########################
sub printOut {
	my ($message, $mlevel) = @_;
	if ($mlevel <= $vlevel){
		print "\n$message\n";
	}
	######################
	sub checkInt {
		my $number = shift;
		if ($number =~ /[^0-9]/){
			return();
		}
		else{
			return($number);
		}
	}
}
###########################
sub helpMessage {
	my $helpmessage = "
	${bold}YOU ARE RUNNING $version on $hostname$norm

This program is freely distributed under a GPL.
Copyright (c) GRL limited: portions of the code are from separate copyrights

\n${bold}USAGE:${norm} oneSeq.pl -sequence_file=<> -seqId=<>  -seqName=<> -refSpec=<> -minDist=<> -maxDist=<> [OPTIONS]

${bold}OPTIONS:$norm

${bold}GENERAL$norm

-h
	Invoke this help method
-version
	Print the program version
-showTaxa
        Print availible Taxa (dependent on the on/off status of database mode)

${bold}REQUIRED$norm

-sequence_file=<>
	Specifies the file containing the seed sequence (protein only) in fasta format. 
	If not provided the program will ask for it.
-seqId=<>
	Specifies the sequence identifier of the seed sequence in the reference protein set. 
	If not provided, the program will attempt to determin it automatically.
-refSpec
	Determines the reference species for the hamstr search. It should be the species the seed sequence was derived from. 
	If not provided, the program will ask for it.
-minDist=<>
	specify the minimum systematic distance of primer taxa for the core set compilation.
	If not provided, the program will ask for it.
-maxDist=<>
	specify the maximum systematic distance of primer taxa to be considered for core set compilation.
	If not provided, the program will ask for it.
-coreOrth=<>
	Specify the number of orthologs added to the core set.

${bold}USING NON-DEFAULT PATHS$norm

-outpath=<>
	Specifies the path for the output directory. Default is $outputPath;
-hmmpath=<>
	Specifies the path for the core ortholog directory. Default is $coreOrthologsPath

${bold}ADDITIONAL OPTIONS$norm

-seqName=<>
        Specifies a name for the search. If not set a random name will be set.
-db
	Run oneSeq.pl in database mode. Requires a mySql database. Only for internatl use.
-filter=[T|F]
	Switch on or off the low complexity filter for the blast search. Default: T
-silent
	Surpress output to the command line
-coreTaxa=<>
	You can provide a list of primer taxa that should exclusively be used for the compilation
	of the core ortholog set
-strict
	Run the final HaMStR search in 'strict mode'. An ortholog is only then accepted when the reciprocity is fulfilled
	for each sequence in the core set.
-force
	Force the final HaMStR search to create output file. Existing files will be overwritten.
-coreStrict
	Run the HaMStR for the compilation of the core set in strict mode.
-checkCoortholgsRef
	During the final HaMStR search, accept an ortholog also when its best hit in the reverse search is not the 
	core ortholog itself, but a co-ortholog of it.
-CorecheckCoorthologsRef
	Invokes the 'checkCoorthologsRef' behavior in the course of the core set compilation.
-rbh
	Requires a reciprocal best hit during the HaMStR search to accept a new ortholog.
-eval_blast=<>
	This option allows to set the e-value cut-off for the Blast search. Default: 1E-5
-eval_hmmer=<>
	This options allows to set the e-value cut-off for the HMM search. Default: 1E-5
-eval_relaxfac=<>
	This options allows to set the factor to relax the e-value cut-off (Blast search and HMM search) for the final HaMStR run. Default: 10
-hitLimit=<>
	Provide an integer specifying the number of hits of the initial pHMM based search that should be evaluated
	via a reverse search. Default: 10
-coreHitLimit=<>
	Provide an integer specifying the number of hits of the initial pHMM based search that should be evaluated
	via a reverse search. Default: 3
-coreOnly
	Set this flag to compile only the core orthologs. These sets can later be used for a stand alone HaMStR search. 
-reuse_core
	Set this flag if the core set for your sequence is already existing. No check currently implemented.
-blast
	Set this flag to determine sequence id and refspec automatically. Note, the chosen sequence id and reference species
	does not necessarily reflect the species the sequence was derived from.
-rep
	Set this flag to obtain only the sequence being most similar to the corresponding sequence in the core set rather
	than all putative co-orthologs.
-coreRep
	Set this flag to invoke the '-rep' behaviour for the core ortholog compilation.
-cpu
	Determine the number of threads to be run in parallel
-batch=<>
	Currently has NO functionality.
-group=<>
	I think it allows to limit the search to a certain systematic group. NOTE: Needs to be checked.
-cleanup
        Temporary output will be deleted.

${bold}SPECIFYING FAS SUPPORT OPTIONS$norm

-fasoff
        Turn OFF FAS support. Default is ON.
-corefilteroff
        Turn OFF the filtering core orthologs. Default is ON.
-minScore=<>
        Specify the threshold for corefilter. Default is 0.8.
-local
        Specify the alignment strategy during core ortholog compilation. Default is local.
-global
        Set the alignment strategy during core ortholog compilation to global.
-countercheck
        Set this flag to counter-check your final profile. The FAS score will be computed in two ways (seed vs. hit and hit vs. seed).

${bold}SPECIFYING EXTENT OF OUTPUT TO SCREEN$norm

-debug
	Set this flag to obtain more detailed information about the programs actions
-silent
	Surpress output to screen as much as possbile

\n\n";
	return($helpmessage);	

}

