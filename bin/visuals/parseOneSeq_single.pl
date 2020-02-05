#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use IO::Handle;
use File::Path;
use Getopt::Long;
use File::Basename;
use Env qw(ONESEQDIR);

#######################
# NAME:         parseOneSep.pl
# AUTHOR:       Vinh Tran, tran@bio.uni-frankfurt.de
# DESCRIPTION:  parsing 1 output of oneSeq.pl (extended.profile) to create input to phyloprofile app
# DATE:         28.11.2018
# SUPPORT:      sge,qsub
# STATUS:       devo
#######################

my $version = 1.0;

sub usage {
    my $msg = shift;
    print "example: perl parseOneSeq.pl -i oneseqOutput.extended.profile -o output.phyloprofile\n";
    print "-i\thamstr oneseq output file (output.extended.profile)\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_o);
getopts('i:o:');

# sanity checks;
my $oneseqFile = ($opt_i) ? $opt_i : usage("ERROR: No oneseq output file given\n");
my $out = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");

### MAIN
unless(-e $oneseqFile){
	usage("ERROR: No extended.profile file found!\n");
}

open(OUT,">$out") || die "Cannot open $out!\n";
print OUT "geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n";

open(IN,$oneseqFile) || die "Cannot open $oneseqFile!\n";
my @in = <IN>;
close (IN);

### get taxon ID, ortho ID and FAS scores from each oneseq output file
foreach my $line(@in){
	chomp($line);	# STRADB|ANOGA@7165@1|Q7Q3C2|1    0.999591652     0.66387534112
	$line =~ s/\|/\t/g;
	my @tmp = split(/\t/,$line);

	my @ncbiID = split(/@/,$tmp[1]);
	my $orthoID = $tmp[0]."|".$tmp[1]."|".$tmp[2]."|".$tmp[3];

	my $fas_F = "0.0";
	if($tmp[4] > 0){$fas_F = $tmp[4];}
	my $fas_B = "0.0";
    if(defined $tmp[5] && $tmp[5] ne "0") {
        $fas_B = $tmp[5];
    }

	print OUT "$tmp[0]\tncbi$ncbiID[1]\t$orthoID\t$fas_F\t$fas_B\n";
}

close (OUT);
print "Finished! Check output file $out\n";
exit;
