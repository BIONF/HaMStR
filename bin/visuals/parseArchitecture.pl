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
# NAME:         parseArchitecture.pl
# AUTHOR:       Vinh Tran, tran@bio.uni-frankfurt.de
# MODIFIED:     Holger Bergmann, bergmann@bio.uni-frankfurt.de
# DESCRIPTION:  parsing output of oneSeq.pl and FAS to create input to phyloprofile app
# DATE:         16.12.2016
# SUPPORT:      sge,qsub
# STATUS:       devo

#######
#SETUP
#######
my $version = 1.0;
my $configure = 1;
my $useIDasis = 1;
if ($configure == 0){
	die "\n\n$version\n\nPLEASE RUN THE CONFIGURE OR CONFIGURE_MAC SCRIPT BEFORE USING parseOneSep.pl\n\n";
}
####################
my $path=$ONESEQDIR;
if (!(defined $path) or !(-e $path)) {
	die "Please set the environmental variabel ONESEQDIR\n";
}  
$path =~ s/\/$//;

## global variables
my $inFile;
my $proFile;
my $groupID;
my $outFile;
my $help;
my $getversion;

##### Command line options
GetOptions ("h"             => \$help,
            "v"             => \$getversion,
            "input=s"       => \$inFile, 
            "profile=s"     => \$proFile,
            "group=s"       => \$groupID,
            "outfile=s"     => \$outFile);


if($help) {
    usage();
    exit;
}
if ($getversion){
    print "You are running $version\n";
    exit;
}

sub usage {
    my $msg = shift;
    print "example: perl parseArchitecture.pl -i xmlOutput.fas -p groupID.extended.profile -g groupID -o output.file\n";
    print "-i\tInput of XML file from FAS calculation\n";
    print "-p\tInput of *extended.profile from hamstr oneSeq\n";
    print "-g\tGroup ID\n";
    print "-o\tOutput file\n";
    die $msg."\n";
}

# global variables
#our($opt_i,$opt_p,$opt_g,$opt_o);
#getopts('i:p:g:o:');

# sanity checks;
#my $inFile = ($opt_i) ? $opt_i : usage("ERROR: No input file given\n");
#my $proFile = ($opt_p) ? $opt_p : usage("ERROR: No input profile given\n");
#my $groupID = ($opt_g) ? $opt_g : usage("ERROR: No group ID given\n");
#my $outFile = ($opt_o) ? $opt_o : usage("ERROR: No output file given\n");

open(OUT0,">".$outFile."_0.domains") || die "Cannot create $outFile!\n";
open(OUT1,">".$outFile."_1.domains") || die "Cannot create $outFile!\n";

open(IN,$inFile) || die "Cannot open $inFile!\n";

### MAIN
my @in = <IN>;
close (IN);
my $in = join("",@in);

### read fitting IDs from profile
open(PRO,$proFile) || die "Cannot open $proFile!\n";
my %taxonmap = ();
while(<PRO>){
	my  $curli = $_;
	chomp($curli);
	my @oneLine = split(/\t/,$curli);
	$taxonmap{$oneLine[0]} = 1; 
}
close(PRO);

### split multiple XML file into individual files if necessary
my @xml = split(/<\/out>/,$in);

foreach my $archi(@xml){
	if(length($archi) > 10){
		my @archiTMP = split(/<\/single_protein>/,$archi); #$archiTMP[0] is seed architecture

		### get search species ID (if available)
		my $searchSpec = "";
		my $direction = "";
		if($archi =~ /(.)+?\.xml/){
			my $hit = $&;		# plaga_4069@5849@1_22941_0_fas.xml
			my @hit = split(/_/,$hit); 
			pop(@hit);	# remove _fas.xml
			$direction = pop(@hit);
			$searchSpec = join("_",@hit);
                        my @stv = split(/@/,$searchSpec);
                        my @v_r = split(/_/,$stv[2]);
                        my $v = $v_r[0];
                        $searchSpec = join("@",($stv[0],$stv[1],$v));
                        #print $searchSpec."\n";
		}

		### get seed ID
		my $seedID = "";
		if($archiTMP[0] =~ /single_protein id=\"(.)+?\"/){
			$seedID = $&;
			$seedID =~ s/single_protein id=//; $seedID =~ s/\"//g;
		}

		### go through all search proteins in this xml file
		my @searchProts = split(/<\/set_protein>/,$archiTMP[1]);

		foreach my $searchProt(@searchProts){
			### get search protein ID(s)
			my $searchID = "";
			if($searchProt =~ /set_protein id=\"(.)+?\"/){
				$searchID = $&;
				$searchID =~ s/set_protein id=//;
                                $searchID =~ s/\"//g;
				if(length($searchSpec) > 0 && !$useIDasis){$searchID = $groupID."|".$searchSpec."|".$searchID;}

				### get info
				my @info = split(/<\/architecture>/,$searchProt);

				### get protein's domain positions
				my $searchDomain; 
				if($direction == 0){
                                        $searchDomain = getDomainPos($groupID,$searchID,$seedID,$info[1],1,$searchSpec);
					print OUT0 $searchDomain;
				} elsif($direction == 1){
                                        $searchDomain = getDomainPos($groupID,$searchID,$seedID,$info[0],1,$searchSpec);
					print OUT1 $searchDomain;
				}
				
				### get seed's best path
				my $seedDomain; 
				if($direction == 0){
                                        $seedDomain = getDomainPos($groupID,$searchID,$seedID,$archiTMP[0],0,$searchSpec);
					print OUT0 $seedDomain;
				} elsif($direction == 1){
                                        $seedDomain = getDomainPos($groupID,$searchID,$seedID,$info[1],0,$searchSpec);
					print OUT1 $seedDomain;
				}
			}
		}
	}
}
close (OUT0);
close (OUT1);

print "Finished! Check output at\n\t",$outFile,"_0.domains\n\t",$outFile,"_1.domains\n";
exit;

sub getDomainPos{
	my ($groupID,$searchID,$seedID,$block,$order,$soi) = @_;

	my @features = split(/feature/,$block);
	my $result = "";

	foreach my $feature (@features){
		if($feature =~ /start/){
			my @info = split(/\n/,$feature);

			my $firstLine = shift(@info);
			my $type = "";
			if($firstLine =~ /type=\"(.)+?\"/){
				$type = $&;
				$type =~ s/type=//; $type =~ s/\"//g;
			}

			my $weight = "NA";
			if($firstLine =~ /weight=\"(.)+?\"/){
				$weight = $&; $weight =~ s/weight=//; $weight =~ s/\"//g;
			}

			foreach my $infoLine(@info){
				chomp($infoLine);
				if($infoLine =~ /start/){
#					print $line,"\n";
					my $start = ""; my $end = "";
					if($infoLine =~ /start=\"\d+\"/){
						$start = $&; $start =~ s/start=//; $start =~ s/\"//g;
					}
					if($infoLine =~ /end=\"\d+\"/){
						$end = $&; $end =~ s/end=//; $end =~ s/\"//g;
					}

                                        ## map IDs
                                        my $exp_ID = "$groupID|$soi|$searchID";
                                        my $rep = $exp_ID."|1";
                                        my $co  = $exp_ID."|0";
                                        my $usedID = "";
                                        if ($taxonmap{$exp_ID}){
                                            $usedID = $exp_ID;
                                        }elsif($taxonmap{$rep}){
                                            $usedID = $rep;
                                        }elsif($taxonmap{$co}){
                                            $usedID = $co;
                                        }else{
                                            usage("ERROR: Given IDs do not match for used ID $usedID !\n")
                                        }

					if($order == 1){
						$result .= "$groupID#$usedID#$seedID\t$usedID\t$type\t$start\t$end\t$weight\n";
					} elsif($order == 0){
						$result .= "$groupID#$usedID#$seedID\t$seedID\t$type\t$start\t$end\t$weight\n";
					}
				}
			}
		}
	}
	return $result;
}
