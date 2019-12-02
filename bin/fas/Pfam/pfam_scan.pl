#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

##################
# This script scans the Pfam-A.hmm database with given input sequences in
# multi-fasta-format. It also looks for further information like clan-ID in the
# Pfam-A.hmm.dat file and creates an output file: pfam_scan.out
#
# HMMR3 (hmmscan), the Pfam-A.hmm database with binary 
# files and also the Pfam-A.hmm.dat file are required.
#
# usage: shell% perl pfam_scan.pl <path to query>
##################
##### Modifications
### added the option to specify the number of CPUs
my $cpu = 4;
########## getting ready ########
my $Directory = cwd();      # save current path/directory
#my $Qquery = shift;       # path to a multi-fasta file with sequences to annotate
my $pfamA_DB = $Directory."/Pfam-hmms";     # path to Pfam-A.hmm
my $outputDir = $Directory."/output_files"; # path to folder where output files are going to be created

my @InputIDs;
my %clans;
########## start main ##########
#main($Qquery);

sub main {
    #my $query = shift;       # path to a multi-fasta file with sequences to annotate
    my ($query,$taxaName) = @_;       # :)
    my $processID = $$;
    ## check everything
    check($query);
    ## scan database
    my $outFile = hmmscan($query,$taxaName,$processID);
    ## get the sequence ID's of the query file
    @InputIDs = getInputIDs($query);
    %clans = getClanIDAndNumber();
    ## analyse the output of hmmscan
    print "Creating output file...\n";
    my $scanresult_fileName = scanAnal($outFile,$taxaName,$processID);

	#### comment the next two lines if you do not want to delete tmp files created by hmmscan
	my $delCommand = "rm -f $outputDir/$outFile";
	system($delCommand);

    print "finished.\n";
	## return output file created by scanAnal()
    return ($scanresult_fileName);
}

########## start subs ##########
sub getClanIDAndNumber {    # write every clan name and the number of members into a hash 
    open (DATFILE, $pfamA_DB."/Pfam-A.hmm.dat")
        or die ("ERROR: could not find file\n");
    my @DATCONTENT = <DATFILE>;
    close DATFILE;

    my %clans;
    for (my $i=32;$i<@DATCONTENT;$i+=8){
        while(!($DATCONTENT[$i]=~/#=GF CL/)){$i++;};
        if($DATCONTENT[$i]=~/#=GF CL/){
            my @clanName = split(/ +/, $DATCONTENT[$i]);
            chomp($clanName[2]);
            if(defined($clans{$clanName[2]})){   # already exists
                my $actualClanNumber = $clans{$clanName[2]};
                $clans{$clanName[2]} = $actualClanNumber+1;
            }
            else{$clans{$clanName[2]} = 1}
        }
    }

    return(%clans);
}

sub scanAnal {
    ### this function analyzes the output of hmmscan extracts all needed informations
    ### and creates "pfam_scan.out"
    my ($outFile,$taxaName,$pid) = ($_[0],$_[1],$_[2]);
    my @extInfos;   # output array content : >headline, seq_id, infos, infos.., >headline, seq_id, infos... 
    my $headline = "## > queryID
## # DomainFamily | type | PfamClan | description
## # hmm_acc|align_start|align_end|env_start|env_end|hmm_start|hmm_end|hmm_len|bit_score|family e-val|instance e-val|numb_clan_memb.|pred_act._site|";
 
    ## open file ready only
    open (SCANOUT, "<".$outputDir."/".$outFile) or die("ERROR: couldnt find $outFile!\n");
    my @ScanOut = <SCANOUT>;
    close SCANOUT;             
    #debugPrint(\@ScanOut);
    my $tempID;     # temporary save for query-prot-ID
    my @tempInfoD;  # temporary save for information of domain
    my @datInfo;    # save space for information of *.dat file
####### loop start
    push(@extInfos, $headline."\n"); # save headline in output array
    for(my $i=11;$i<@ScanOut;$i++){ # reading in the lines of hmmscan.out and getting out information
        if($ScanOut[$i]=~/Query:/){ # look for line that starts with "Query:"
            my @getID = split(/ +/, $ScanOut[$i]);
            my %fam_seq_info;
            #print "ID: ".$getID[1]."\n"; # DEBUG
            #print "1. ".$ScanOut[$i]."\n";  # DEBUG
            $getID[1]="$getID[1]";
            my $y=$i;
            my $noHits=0;
            while($y<$i+10){
                if($ScanOut[$y]=~/\[No hits detected that satisfy reporting thresholds\]/){
                        $noHits = 1;
                        $tempID = shift(@InputIDs);
                        $tempID =~ s/\r//;
                        #print "1 x".$getID[1]."x eq x".$tempID."x \n";  # DEBUG
                        if($getID[1] eq $tempID){    # check if Query is InputID
                        #print "works! 1 \n"; # DEBUG
                            push(@extInfos, "\n>".$tempID."\n");  # save InputID in string
                            push(@extInfos, "# [No hits detected that satisfy reporting thresholds\]\n");
                        }else{print "problem 1!\n"} # DEBUG
                        last;
                }$y++;
            }
            if($noHits==0){
                    $tempID = shift(@InputIDs);
                    $tempID =~ s/\r//;
                    #print "1 x".$getID[1]."x eq x".$tempID."x \n";  # DEBUG
                    if($getID[1] eq $tempID){    # check if Query is InputID
                        #print "works! 2 \n"; # DEBUG
                        push(@extInfos, "\n>".$tempID."\n");  # save InputID in string
                        #$i+=8;
                    }else{print "problem 2!\n"} # DEBUG

                    ## reading complete sequence information
                    while(!($ScanOut[$i]=~/Domain annotation for each model:/)){
                        if ($ScanOut[$i]=~/E-value  score  bias/){
                            $i = $i+2;
                            while($ScanOut[$i] ne "\n"){#proceed until empty line in ScanOut
                                my @c_seq_info = split(/ +/, $ScanOut[$i]);
                                $fam_seq_info{$c_seq_info[9]} = $c_seq_info[1];
                                $i++;
                            }
                        }                        
                        $i++;
                    }

                    ## reading domain annotation for each model (continue with instancewise information)
                    if($ScanOut[$i]=~/Domain annotation for each model:/){  # if this is domain annot. line 
                        my $whileExit=0;    # while this variable is "0" the while loop will be repeated
                        while($whileExit==0){   
                            # loop repeats getting "domain annot." info. and also infos out of the *.dat file
                            # the infos are saved into strings and then pushed to the output array.
                            # Each found hit for one query will be analyzed in this loop.
                            my $infoString;
                            $i++;   # go to next line
                            @tempInfoD = split(/ +/, $ScanOut[$i]);  # split this line (>> hmm name description)

                            # datInfo = AC(hmm_acc), ID(hmm_name), TP(type), ML(length), CL(clan)
                            @datInfo = getDatFileInfo($tempInfoD[1]);
                            if(!defined($datInfo[0])){$datInfo[0]="---"};
                            if(!defined($datInfo[1])){$datInfo[1]="---"};
                            if(!defined($datInfo[2])){$datInfo[2]="---"};
                            if(!defined($datInfo[3])){$datInfo[3]="---"};
                            if(!defined($datInfo[4])){$datInfo[4]="---"};
                            my $tempHeader; # variable to save header information
                            $tempHeader.= "# family: ".$datInfo[1]." | type: ".$datInfo[2]." | clan: ".$datInfo[4]." | description: ";

                            for(my $i2=2;$i2<@tempInfoD;$i2++){     # save description-part
                                chomp($tempInfoD[$i2]);
                                $tempHeader .= $tempInfoD[$i2]." ";
                            }
                            $tempHeader.= "| E-value: ".$fam_seq_info{$datInfo[1]} ;

                            unless($ScanOut[$i+1]=~ /\[No individual domains that satisfy reporting thresholds \(although complete target did\)\]/){                    
                                    push(@extInfos, $tempHeader."\n");   # save header in output array


                                    $i+=3;  # jump to line with information of interest#my $escape=0;

                                    while(defined($ScanOut[$i]) && length($ScanOut[$i])>1){# && $ScanOut[$i]!~/No individual domains that satisfy reporting thresholds (although complete target did)/){
                                        # @tempInfoD = #,score,bias,c-Evalue,i-Evalue,hmmfrom,hmm to,alifrom,ali to,envfrom,env to,acc
                                        @tempInfoD = split(/ +/, $ScanOut[$i]);
                                        # now writing all infos into a single string, in the right order
                                        #$tempInfoD[6]  contains domain-/instancewise E-Values
                                        #               contains 
                                        #debugPrint(\@tempInfoD);
                                        $infoString .= "#".$datInfo[0]."|".$tempInfoD[10]."|".$tempInfoD[11]."|".$tempInfoD[13]."|".$tempInfoD[14];
                                        $infoString .= "|".$tempInfoD[7]."|".$tempInfoD[8];
                                        $infoString .= "|".$datInfo[3];
                                        $infoString .= "|".$tempInfoD[3]."|".$fam_seq_info{$datInfo[1]}."|".$tempInfoD[5];
                                        if(defined($clans{$datInfo[4]})){
                                            $infoString .= "|".$clans{$datInfo[4]}."|";
                                        }else{$infoString .= "|--"."|"}

                                        push(@extInfos, $infoString."\n");   # save infos in output array
                                        $infoString="";
                                        $i++;
                                    }
                            }else {
                                #print $ScanOut[$i+1],"\n";
                                $i += 2;
                            }
                            if(!($ScanOut[$i+1]=~/\>\>/)){$i+=10;$whileExit=1;} # now we can jump 16 lines
                        }
                    }
                }
            }
        }
###### loop end

	my $outPut_fileName="pfam_scan_".$taxaName;

	open(OUTPUT, ">$outputDir/$outPut_fileName$pid.out") or die("ERROR: could not write output file.\n");
	foreach(@extInfos){
		print OUTPUT $_;
	}
	close OUTPUT;

	return("$outPut_fileName$pid.out");
}

sub debugPrint{
    my @currentarray = @{$_[0]};
    print "\n";
    foreach(@currentarray){
    print $_ . "|";
    }
    print "\n";
}

sub getDatFileInfo {
    ### with the query ID as parameter this function gets out information
    ### of the Pfam-A.hmm.dat file and returns them as array

    my $domainName = shift(@_);
    my @return;
    my @tempdatContent;

    open (DATFILE, $pfamA_DB."/Pfam-A.hmm.dat")   # open file
        or die("ERROR: couldnt find Pfam-A.hmm.dat!\n");
    my @datContent = <DATFILE>;    # save content in array
    close DATFILE;              # close file

    for(my $i=1;$i<@datContent;$i+=8){ # search datFile
        if(!($datContent[$i]=~/\#=GF ID/)){  # if this line NOT begins with "#=GF ID"...
            while(!($datContent[$i]=~/\#=GF ID/)){   # ...increase $i until it does
                $i++;
            }
        }
        if($datContent[$i]=~/$domainName/){ # else if this line contains the searched domainName

            @tempdatContent = split(/ +/, $datContent[$i+1]);  # split the line
            chomp($tempdatContent[2]);
            push(@return, $tempdatContent[2]);  # save first return value -> AC
 
            push(@return, $domainName);  # save second return value -> ID 

            while(!($datContent[$i]=~/\/\//)){
                $i++;
                
                if($datContent[$i]=~/\#=GF TP/){
                    @tempdatContent = split(/ +/, $datContent[$i]);  # split the line
                    chomp($tempdatContent[2]);
                    push(@return, $tempdatContent[2]);  # save third return value -> TP
                }
                if($datContent[$i]=~/\#=GF ML/){
                    @tempdatContent = split(/ +/, $datContent[$i]);  # split the line
                    chomp($tempdatContent[2]);
                    push(@return, $tempdatContent[2]);  # save fourth return value -> ML
                }
                if($datContent[$i]=~/\#=GF CL/){
                    @tempdatContent = split(/ +/, $datContent[$i]);  # split the line
                    chomp($tempdatContent[2]);
                    push(@return, $tempdatContent[2]);  # save last return value -> CL
                }
            }
            last;
        }
    }
    # return = AC(hmm_acc), ID(hmm_name), TP(type), ML(length), CL(clan)
    return(@return);
}

sub getInputIDs {
    my @InputIDs;
    my $query = shift;

    open(QUERY, $query)
        or die ("ERROR: couldnt find query file: $query\n");
    my @query = <QUERY>;
    close QUERY;

    foreach(@query){
        if($_=~ /\>/){
            $_=~ s/\>//;
            chomp($_);
            push(@InputIDs, $_)
        }
    }

    return(@InputIDs);
}

sub check {
    my $query = shift;
    ### this function checks if all requirements are fullfilled
    if(!(-d ($pfamA_DB))){
        mkdir($pfamA_DB);
    };

    if(!(-d ($outputDir))){
        mkdir($outputDir);
    };
    ## check if all required files exist
    # check if database exists
    if(!(-e ("$pfamA_DB/Pfam-A.hmm"))){
        die("ERROR: Database \"Pfam-A.hmm\" is missing, please download and/or put into folder: $pfamA_DB \n");
    };
    # check if .dat file exists
    if(!(-e ("$pfamA_DB/Pfam-A.hmm.dat"))){
        die("ERROR: Pfam-A.hmm.dat is missing, please download and/or put into folder: $pfamA_DB \n");
    };
    # check if binary files exist
    if(!(-e ("$pfamA_DB/Pfam-A.hmm.h3f") and -e ("$pfamA_DB/Pfam-A.hmm.h3i") and -e ("$pfamA_DB/Pfam-A.hmm.h3m") and -e ("$pfamA_DB/Pfam-A.hmm.h3p"))){
        print "ERROR: Binary files: Pfam-A.hmm.h3f, *.h3i, *.h3m, *.h3m are missing! Create now? [Y/n] ";
        my $choose = <STDIN>;
        chomp($choose);

        if(($choose eq "y") or ($choose eq "Y")){
            print "Creating binary files...\n";
            system("hmmpress -f $pfamA_DB/Pfam-A.hmm");
            print "Finished creating binary files.\n";
        }
        else{die("ERROR: Binary files: Pfam-A.hmm.h3f, *.h3i, *.h3m, *.h3m are missing!\n");};
    };
    ## check if input file exists
    if(!defined($query) or !(-e ($query))){
        die("ERROR: Input file does not exist or could not be found: $query\n");
    };
}

sub hmmscan {
### this function uses hmmscan
    #my $query = shift;
    my ($query,$taxaName,$pid) = ($_[0],$_[1],$_[2]);
    #my $outPut_fileName = "hmmscan";
    my $outPut_fileName = "hmmscan_".$taxaName;
    print "Starting hmmscan now...\n";
    ## in the next line hmmscan is used to find domains in the pfam database.
    ## the output file is in pfam_format
    # TODO:	--domtblout <f>
    # VAL:	--E <E-value> 		(full sequence length)
    # VAL:	--domE <c-Evalue>	(domainwise/instancewise)
    system("hmmscan --cpu $cpu -E 0.01 --domE 0.1 --noali -o \"$outputDir/$outPut_fileName$pid.out\" \"$pfamA_DB\"/Pfam-A.hmm \"$query\""); #--pfamtblout
    print "hmmscan finished.\n";

    return($outPut_fileName.$pid.".out");
}

1;
