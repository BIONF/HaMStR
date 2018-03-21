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
# AUTHOR:		Holger Bergmann, bergmann@bio.uni-frankfurt.de
# DESCRIPTION:  parsing output of oneSeq.pl and FAS to create input to phyloprofile app
# DATE:         16.12.2016

#######
#SETUP
#######
my $version = 1.3;
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
my $debug = 0;

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
    print "example: perl parseArchitecture.pl -i scores_1_fas.collection -p groupID.extended.profile -g groupID -o output.file\n";
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

my ($in_base, $in_path, $in_suffix) = fileparse( $inFile, qr/\.[^.]*/ );
my @direction = split(/_/,$in_base);
#print $direction[1]."\n";


open(OUT,">".$outFile."_".$direction[1].".domains") || die "Cannot create $outFile!\n";

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
### split by endnote of <architectures> section
my @xml = split(/<\/architectures>/,$in);
if ($debug){print $xml[1]."\n";}
#---------------> further parsing in abhaengigkeit von direction
#---------------> 1: template = seed, query = ortholog
#---------------> 0: template = ortholog, query = seed

foreach my $archi(@xml){
    
    if(length($archi) > 10){
        if ($direction[1] == "1"){
            my @archiTMP = split(/<\/template>/,$archi); 

            ### get information #
            my $seed        = "";
            my $query       = "";
            my $seedlen		= "";
            my $querylen	= "";
            my $queryid     = "";
            my $queryID     = "";
            my $direction   = "";
            if($archi =~ /(.)+?\.xml/){
                my $hit = $&;
                my @hit = split(/_/,$hit);
                pop(@hit);
                $direction = pop(@hit); #better save than sorry
                if ($direction[1] != $direction){
                    print "\nWARNING: please check file format for:\t".$inFile."\n";
                }
           }
           if($archiTMP[0] =~ /query id=\"(.)+?length=\"(.)+?\"/){
                my $hit = $&;
                $hit =~ s/query id=//; $hit =~ s/\slength=.*//; $hit =~ s/\"//g;
                my @hit =  split(/\|/,$hit);
                $queryid = $hit[1];

                $query = $hit[0];
                if ($debug){print $query." and ".$queryid."\ndirection: ".$direction."\n";}
            }

            if($archiTMP[0] =~ /template id=\"(.)+?length=\"(.)+?\"/){
            	my $templateline = $&;
            	$seed = $templateline;
            	$seed =~ s/template id=//; $seed =~ s/\sscore=.*//; $seed =~ s/\"//g;
                if ($debug){print "SEED ID: ".$seed."\n";}
                
               	$seedlen = $templateline;
                $seedlen =~ s/.*length=//;
   	            $seedlen =~ s/\"//g;
   	            if ($debug){print "SEED LENGTH: ".$seedlen."\n";}
            }
            if($archiTMP[0] =~ /query id=\"(.)+?length=\"(.)+?\"/){
            	my $queryline = $&;
               	$querylen = $queryline;
                $querylen =~ s/.*length=//;
   	            $querylen =~ s/\"//g;
   	            if ($debug){print "QUERY LENGTH: ".$querylen."\n";}
            }
            
            if(length($query) > 0 && !$useIDasis){
                $queryID = $groupID."|".$query."|".$queryid;
            }else{
                $queryID = $queryid;
            }
            if ($debug){print "queryID: ".$queryID."\n";}

            ## information in fields of @archiTMP
            #$archiTMP[0]: paths (seed and query)
            #$archiTMP[1]: seed architecture
            #$archiTMP[2]: query architecture

            ### get seed path
            #$fixture[0]: seed path
            #$fixture[1]: query path
            my @fixture = split(/<\/template_path>/,$archiTMP[0]);
            my %seedPath = getPathInfo($groupID,$queryID,$seed,$fixture[0],0,$query);
            my %queryPath = getPathInfo($groupID,$queryID,$seed,$fixture[1],1,$query);
            
#            print "seed\n";
#            foreach my $key(keys %seedPath){
#                print $key." --> ".$seedPath{$key}."\n";                
#            }
#            print "query\n";
#            foreach my $yek(keys %queryPath){
#                print $yek." --> ".$queryPath{$yek}."\n";                
#            }

            ### get query domains + path
            my $queryDomains = getDomainPos($groupID,$queryID,$seed,$archiTMP[2],1,$query, $querylen, \%queryPath);
            print OUT $queryDomains;
            
            ### get seed domains + path
            my $seedDomains = getDomainPos($groupID,$queryID,$seed,$archiTMP[1],0,$query, $seedlen, \%seedPath);
            print OUT $seedDomains;
        }
        if ($direction[1] == "0"){
            my @archiTMP = split(/<\/template>/,$archi); 
            ### get information #
            my $seed        = "";
            my $query       = "";
            my $seedlen		= "";
            my $querylen	= "";
            my $queryid     = "";
            my $queryID     = "";
            my $direction   = "";
            if($archi =~ /(.)+?\.xml/){
                my $hit = $&;
                my @hit = split(/_/,$hit);
                pop(@hit);
                $direction = pop(@hit); #better save than sorry
                if ($direction[1] != $direction){
                    print "\nWARNING: please check file format for:\t".$inFile."\n";
                }
            }
            if($archiTMP[0] =~ /template id=\"(.)+?length=\"(.)+?\"/){
                my $hit = $&;
                $hit =~ s/template id=//; $hit =~ s/\sscore=.*//; $hit =~ s/\"//g;
                my @hit =  split(/\|/,$hit);
                $queryid = $hit[1];

                $query = $hit[0];
                if ($debug){print $query." and ".$queryid."\ndirection: ".$direction."\n";}
            }   

            if($archiTMP[0] =~ /query id=\"(.)+?length=\"(.)+?\"/){
            	my $queryline = $&;
            	$seed = $queryline;
            	$seed =~ s/query id=//; $seed =~ s/\slength=.*//; $seed =~ s/\"//g;
                if ($debug){print "SEED ID: ".$seed."\n";}
                
               	$seedlen = $queryline;
                $seedlen =~ s/.*length=//;
   	            $seedlen =~ s/\"//g;
   	            if ($debug){print "SEED LENGTH: ".$seedlen."\n";}
            }
            if($archiTMP[0] =~ /template id=\"(.)+?length=\"(.)+?\"/){
            	my $templateline = $&;
               	$querylen = $templateline;
                $querylen =~ s/.*length=//;
   	            $querylen =~ s/\"//g;
   	            if ($debug){print "QUERY LENGTH: ".$querylen."\n";}
            }

            if(length($query) > 0 && !$useIDasis){
                $queryID = $groupID."|".$query."|".$queryid;
            }else{
                $queryID = $queryid;
            }
            if ($debug){print "queryID: ".$queryID."\n";}

            ## information in fields of @archiTMP
            #$archiTMP[0]: paths (query and seed)
            #$archiTMP[1]: query architecture
            #$archiTMP[2]: seed architecture

            ### get seed path
            #$fixture[0]: query path
            #$fixture[1]: seed path
            # please note, that seed&query informations are labeled differently in case "0"
            my @fixture = split(/<\/template_path>/,$archiTMP[0]);
            my %seedPath = getPathInfo($groupID,$queryID,$seed,$fixture[1],0,$query);
            my %queryPath = getPathInfo($groupID,$queryID,$seed,$fixture[0],1,$query);

            ### get query domains
            my $queryDomains = getDomainPos($groupID,$queryID,$seed,$archiTMP[1],1,$query, $querylen,\%queryPath);
            print OUT $queryDomains;

            my $seedDomains = getDomainPos($groupID,$queryID,$seed,$archiTMP[2],0,$query, $seedlen,\%seedPath);
            print OUT $seedDomains;
        }
    }
}
close (OUT);

print "Finished! Check output at\n\t",$outFile,"_".$direction[1].".domains\n";
exit;

sub getDomainPos{
    my ($groupID,$searchID,$seedID,$block,$order,$soi,$seqlen,%pathinfo) = ($_[0], $_[1], $_[2], $_[3], $_[4], $_[5], $_[6], %{$_[7]});

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
            foreach my $infoLine(@info){
                chomp($infoLine);
                if($infoLine =~ /start/){
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
                    ### print query infos
                    if($order == 1){
                        my $featureinfo = $type.$start.$end;
                        if (exists $pathinfo{$featureinfo}){
                            $result .= "$groupID#$usedID\t$usedID\t$seqlen\t$type\t$start\t$end\t$pathinfo{$featureinfo}\tY\n";
                        }else{
                            $result .= "$groupID#$usedID\t$usedID\t$seqlen\t$type\t$start\t$end\t$weight\tN\n";
                        }
                    #print seed infos        
                    } elsif($order == 0){
                        my $featureinfo = $type.$start.$end;
                        if (exists $pathinfo{$featureinfo}){
                            $result .= "$groupID#$usedID\t$seedID\t$seqlen\t$type\t$start\t$end\t$pathinfo{$featureinfo}\tY\n";
                        }else{
                            $result .= "$groupID#$usedID\t$seedID\t$seqlen\t$type\t$start\t$end\t$weight\tN\n";
                        }
                    }
                }
            }
        }
    }
    return $result;
}
sub getPathInfo{
    my ($groupID,$searchID,$seedID,$block,$order,$soi) = @_;

    my @features = split(/feature/,$block);
    my %pathinfo;

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
            if($firstLine =~ /corrected_weight=\"(.)+?\"/){
                $weight = $&; $weight =~ s/corrected_weight=//; $weight =~ s/\"//g;
            }

            foreach my $infoLine(@info){
                chomp($infoLine);
                if($infoLine =~ /start/){
                    if($debug){print $infoLine,"\n";}
                    my $start = ""; my $end = "";
                    if($infoLine =~ /start=\"\d+\"/){
                        $start = $&; $start =~ s/start=//; $start =~ s/\"//g;
                    }
                    if($infoLine =~ /end=\"\d+\"/){
                        $end = $&; $end =~ s/end=//; $end =~ s/\"//g;
                    }
                    my $featureID = $type.$start.$end;
                    $pathinfo{$featureID} = $weight;
                }
            }
        }
    }
    return %pathinfo;
}
