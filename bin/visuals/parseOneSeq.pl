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
# DESCRIPTION:  parsing output of oneSeq.pl to create input to phyloprofile app
# DATE:         05.02.2018
# SUPPORT:      sge,qsub
# STATUS:       devo

#######
#SETUP
#######
my $version = 1.0;
my $configure = 1;
if ($configure == 0){
	die "\n\n$version\n\nPLEASE RUN THE CONFIGURE OR CONFIGURE_MAC SCRIPT BEFORE USING parseOneSep.pl\n\n";
}
####################
my $path=$ONESEQDIR;
if (!(defined $path) or !(-e $path)) {
	die "Please set the environmental variabel ONESEQDIR\n";
}  
$path =~ s/\/$//;


sub usage {
    my $msg = shift;
    print "example: perl parseOneSeq.pl -i oneseqOutput -o output.phyloprofile\n";
    print "-i\thamstr oneseq output folder (contains *.extended.profile)\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_o);
getopts('i:o:');

# sanity checks;
my $oneseqFolder = ($opt_i) ? $opt_i : usage("ERROR: No oneseq output folder given\n");
my $out = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");

### MAIN
my @allOutFiles = glob("$oneseqFolder/*.extended.profile");
unless(@allOutFiles){
	usage("ERROR: No extended.profile file found!\n");
}

open(OUT,">$out") || die "Cannot open $out!\n";
print OUT "geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n";
foreach my $file(@allOutFiles){
	open(IN,$file) || die "Cannot open $file!\n";
	my @in = <IN>;
	close (IN);

	### get taxon ID, ortho ID and FAS scores from each oneseq output file
	foreach my $line(@in){
		chomp($line);	# STRADB|ANOGA@7165@1|Q7Q3C2|1    0.999591652     0.66387534112
#		print $line,"\n";
		$line =~ s/\|/\t/g;
		my @tmp = split(/\t/,$line);

		my @ncbiID = split(/@/,$tmp[1]);
		my $orthoID = $tmp[0]."|".$tmp[1]."|".$tmp[2]."|".$tmp[3];

		my $fas_F = "0.0";
		if($tmp[4] > 0){$fas_F = $tmp[4];}
		my $fas_B = "0.0";
    if(defined $tmp[5] && $tmp[5] ne "0"){$fas_B = $tmp[5];}

		print OUT "$tmp[0]\tncbi$ncbiID[1]\t$orthoID\t$fas_F\t$fas_B\n";
	}
}

close (OUT);
print "Finished! Check output file $out\n";
exit;
