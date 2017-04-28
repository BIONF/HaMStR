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
# MODIFIED:     Holger Bergmann, bergmann@bio.uni-frankfurt.de
# DESCRIPTION:  parsing output of oneSeq.pl to create input to phyloprofile app
# DATE:         16.12.2016
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
    print "example: perl parseOneSeq.pl -i oneseqOutFolder -o output.matrix\n";
    print "-i\thamstr oneseq output file (*.extended.profile)\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_t,$opt_o);
getopts('i:t:o:');

# sanity checks;
my $oneseqFile = ($opt_i) ? $opt_i : usage("ERROR: No oneseq output file given\n");
my $out = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");

### MAIN
my @allOutFiles = glob("$oneseqFile");
unless(@allOutFiles){
	usage("ERROR: No extended.profile file found!\n");
}

my %taxaList;	# list of all taxa
my %allGenes;	# list of all genes
my %fas0;	# $fas0{$protID#$taxonID} = FAS_score
my %fas1;	# $fas1{$protID#$taxonID} = FAS_score

foreach my $file(@allOutFiles){
	open(IN,$file) || die "Cannot open $file!\n";
	my @in = <IN>;
	close (IN);

	### get taxon ID, ortho ID and FAS scores from each oneseq output file
	foreach my $line(@in){
		chomp($line);	# Arath01153|aquco_5436@218851@1|28258|1	0.99512260472
#		print $line,"\n";
		$line =~ s/\|/\t/g;		# Arath01153	aquco_5436@218851@1	28258	1	0.99512260472
		my @tmp = split(/\t/,$line);

		my $geneID = $tmp[0];
		my $fas0 = "NA";
		my $fas1 = "NA";
		if($tmp[4]){ $fas1 = $tmp[4];}
		if($tmp[5]){ $fas0 = $tmp[5];}
		
		my @hit = split(/\@/,$tmp[1]);
		my $taxonID = $hit[1];
#		my $hitID = $hit[0].":".$tmp[2];
		my $hitID = $tmp[0]."|".$tmp[1]."|".$tmp[2]."|".$tmp[3];
#		print "$geneID - $taxonID - $fas";<>;

		### save to %taxaList, %allGenes and %fas
		$taxaList{"ncbi$taxonID"} = 1;
		unless($fas0{"$geneID#ncbi$taxonID"}){
			$fas0{"$geneID#ncbi$taxonID"} = $hitID."#".$fas0;
			$fas1{"$geneID#ncbi$taxonID"} = $hitID."#".$fas1;
		} else {
			if($fas0{"$geneID#ncbi$taxonID"} =~ /\#NA/){
				$fas0{"$geneID#ncbi$taxonID"} = $hitID."#".$fas0;
			} else {
				unless($fas0 eq "NA"){
					my @fasTMP = split(/\#/,$fas0{"$geneID#ncbi$taxonID"});		
					if($fasTMP[1] < $fas0){
						$fas0{"$geneID#ncbi$taxonID"} = $hitID."#".$fas0;
					}
				}
			}
			
			if($fas1{"$geneID#ncbi$taxonID"} =~ /\#NA/){
				$fas1{"$geneID#ncbi$taxonID"} = $hitID."#".$fas1;
			} else {
				unless($fas1 eq "NA"){
					my @fasTMP = split(/\#/,$fas1{"$geneID#ncbi$taxonID"});		
					if($fasTMP[1] < $fas1){
						$fas1{"$geneID#ncbi$taxonID"} = $hitID."#".$fas1;
					}
				}
			}
		}
		$allGenes{$geneID} = 1;
	}
}

### print output
open(OUT0,">".$out."_0.matrix") || die "Cannot create $out\_0.matrix!\n";
open(OUT1,">".$out."_1.matrix") || die "Cannot create $out\_1.matrix!\n";

my @allTaxa = sort keys %taxaList;
my $allTaxa = join("\t",@allTaxa);
print OUT0 "geneID\t$allTaxa\n";
print OUT1 "geneID\t$allTaxa\n";

foreach my $gene(sort keys %allGenes){
	print OUT0 $gene;
	print OUT1 $gene;
	
	foreach my $taxon(sort @allTaxa){
		if($fas0{"$gene#$taxon"}){
			print OUT0 "\t",$fas0{"$gene#$taxon"};
		} else {
			print OUT0 "\t","NA";
		}
		
		if($fas1{"$gene#$taxon"}){
			print OUT1 "\t",$fas1{"$gene#$taxon"};
		} else {
			print OUT1 "\t","NA";
		}
	}
	print OUT0 "\n";
	print OUT1 "\n";
}
close (OUT0);
close (OUT1);

print "Finished! Check output file\n\t",$out,"_0.matrix\n\t",$out,"_1.matrix\n";

exit;