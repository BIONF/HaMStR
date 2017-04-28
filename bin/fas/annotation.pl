#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use Getopt::Long;
use File::Path;
use File::Basename;


########## dev. version v6 ############
# This script is using a fasta file as input for different protein annotation 
# tools. Only needed information will then be filtered out of the results of
# every tool and will then be written into an output file as an xml file.
# 
# If a foldername is given a folder will be created in "./annot/foldername" so the
# created annotation xml files can be saved permanently. If no foldername is given the
# created xml files will either be saved in the "./annot/query" folder or the 
# "./annot/proteom" folder, depending on the second command "q or p" and existing
# xml files in these folders will be overwritten.
# 
# usage: shell$ perl annotations.pl TODO
#
# note : either "q or p" or a foldername need to be set. if both are given the second
# argument will be ignored. 
#
# MODIFICATIONS made by holger bergmann:
#	25.05.2015	- missing clan annotations will be filled with the domain family name it self
#				- hereby the clan score will not be 0 if two identical domains (annotated in q and p) have no clan annotated
#	01.06.2015	- added var $fasPath: location of the FACT folder. Default: my $fasPath = "/home/holger/src/FACT/" (must be adapted to your environment)
#	02.06.2015	- changed output format
#				- added sequence length calculation (see sub: seq_id_len()) 
#				- compare ids: check for clean labels
#				- clan info is assigned directly to feature type
#				- example: 
#
#				-	<?xml version="1.0"?>
#					<tool name="PfamDB">
#					        <protein id="three" length="1000">
#					                <feature type="Pfam1" instance="1" clan="Pfam1">
#					                        <region start="2" end="119"/>
#					                </feature>
#					                <feature type="Pfam2" instance="2" clan="Pfam2">
#					                        <region start="140" end="419"/>
#					                        <region start="600" end="825"/>
#					                </feature>
#					        </protein>
#					</tool>

#	Version 4 (06.07.2015)
#			- Creating separate temporary files in SMART and PFAM. Deleting them after use.
#			- suitable for greedyFACT.py
#	Version 5 (modifications suggested by Vinh Tran)
#			- updated smart_scan.pl: minor bug fix
#			- CAST segmentation fault: minor bug fix
#			- updated model for TMHMM2.0: extended wildcard alphabet
#   Version 6 (minor modifications by Holger Bergmann)
#           - tool/database names set up in SETUP TOOL BOX
#           - path and location changed to independent design/nameing conventions
#           - option -force: force redo of annotations
#           - option -extraxt: feature/annotaion extraction of existing annotations (no redo)
#   Version 7 (minor modifications by Holger Bergmann, May 2016)
#           - additional field in output XML (hmm based search): evalue: contains the evalue (family, pfam) of features
#   Version 8 (fixing smart over/under predictions by Holger Bergmann, October 2016)
#           - SMART annotations via hmmscan rather hmmsearch
##########

#### SETUP TOOL BOX ####
my $CAST_tool       = "cast";
my $TMHMM_tool      = "decodeanhmm";
my $COILS_tool      = "COILS2";
my $SIGNALP_tool    = "signalp.pl";
my $SEG_tool        = "seg";
my $PFAM_tool       = "pfam_scan.pl";
my $SMART_tool      = "smart_scan_v4.pl";
## number of above stated programs for annotation
my $tool_count      = 7;
my $version         = 1.0;

#### SETUP PATH ####
my $location	= abs_path($0);
my ($base, $fasPath, $suffix) = fileparse( $location, qr/\.[^.]*/ );
chdir($fasPath);

#### PRINT SETUP ####
print "\n#########################\n";
print "--> Running ".$base.$suffix."\n\n";

#### CHECK OPTIONS ####
my $usage           =	"";
my $helpmessage     = getHelpMessage();
 
my $verbose = 0;
my $help;
my $fasta;
my $annot;
my $qORp;
my $redo = 0;
my $force = 0;      #default
my $extract = '';   #default
my $empty = 0;      #default
my $regular = 0;    #default
# command line arguments:
GetOptions( "h"         => \$help,
            "fasta=s"	=> \$fasta,
            "path=s"	=> \$annot,
            "name=s"    => \$qORp,
            "extract=s" => \$extract,
            "redo=s"    => \$redo,
            "force"     => \$force
);

# help
if (defined $help){
    print $helpmessage;
    print $usage;
    die ( "\n ...Exiting ".$base.$suffix."\n");
}
## check for annotation files for extraction
## if $extract is given not annotations will be made
if ($extract ne ''){
    my @fileList = glob($annot."/*.xml");
    if(scalar(@fileList) == $tool_count){
        ## using annotations and extract given gene_id
        print "Found annotation files:\n\t $annot\n";
        print "--> Extracting annotations for given gene id:\n\t $qORp.\n";
        extractAnnotation($qORp);
        print "\n--> Annotations extracted an stored at:\n\t $extract/\n";
        print "--> Exiting script ".$base.$suffix."\n";
        exit;
    }else{
        print "ERROR: Given annotations do not exist or may be incomplete.\n";
        exit;
    }
}

# create given output directory
my $dirOut = $annot."/";
if(!(-d($dirOut))){
	mkdir($dirOut);
}

if($qORp eq "q"){$dirOut.="query/"}
elsif($qORp eq "p"){$dirOut.="proteom/"}
elsif(length($qORp)>0){$dirOut.="$qORp/"}
else{print("ERROR: choose (q)uery or (p)roteom or set a foldername to save the xml files.\n")}

## check for given output if it already exists (specified by -path option)
## if it exists, check if XMl files exist and skip annoation if they do and option force is not set
if(!(-d($dirOut))){     
    mkdir($dirOut);
    $regular = 1;
        
} else {
    my @fileList = glob($dirOut."/*.xml");
    if(scalar(@fileList) == $tool_count && !$force && !$redo){
	print "Found annotation files:\n\t $dirOut. Using these.\n";
        ## exit statement:
	die("Skipping annotations. Using existing files.\n");
    }elsif (!$force && !$redo){
	print "Directory:\n\t $dirOut already exists, but seems to be empty or incomplete. Starting annotation.\n";
        ## set $empty for later evaluation
        $empty = 1;
    }elsif ($force && !$redo){
	print "Directory:\n\t $dirOut already exists but may be empty or incomplete, annotations will be redone due to option -force.\n";
	## delete existing files in output dir and redo the annotations!
	my $delCommand = "rm -f $dirOut/*";
	system($delCommand);
    }elsif ($redo){
        if ($redo eq "cast" or $redo eq "tmhmm" or $redo eq "coils" or $redo eq "signalP" or $redo eq "seg" or $redo eq "pfam" or $redo eq "smart"){
            print "Directory:\n\t$dirOut already exists,\n\t$redo annotations will be redone due to option -redo=$redo.\n";
            my $delCommand = "rm -f $dirOut/$redo.xml";
            system($delCommand);
        }else {
            die("\nPlease check options. -redo=$redo is not a valid option. Something went wrong. Exiting...\n\n");
        }
        
    }else{
        die("\nPlease check options and output directory. Something went wrong. Exiting...\n\n");
    }
}

#### INIT ####
my %any_ids;	# hash with sequence lengths

#### STOPPING TIME START ####
my @time    = localtime(time());
my $milsec1 = time;
$milsec1    -= int($milsec1);
$milsec1    *= 100;
$milsec1    = int($milsec1);

if($time[1]<10){$time[1] = "0".$time[1]}
$time[4] += 1;
if($time[4]<10){$time[4] = "0".$time[4]}
my $START = $time[3]."/".$time[4]."/".($time[5]+1900)." ".$time[2].":".$time[1].":".$time[0].":".$milsec1;

########## START MAIN ##########
#acquire sequence length and check identifier
print "--> acquiring sequence lengths\n";
seq_id_len($fasta);
print "...done: sequence lengths are calculated for input file.\n";

if (($redo eq 'cast') or $force or $empty or $regular){
    print "--> starting: $CAST_tool\n";
    CASTing();
}

if (($redo eq 'tmhmm') or $force or $empty or $regular){
    print "--> starting: $TMHMM_tool\n";
    TMHMM20();
}

if(($redo eq 'coils') or $force or $empty or $regular){
    print "--> starting: $COILS_tool\n";
    coils();
}

if (($redo eq 'signalp') or $force or $empty or $regular){
    print "--> starting: $SIGNALP_tool\n";
    signalp();
}

if (($redo eq 'seg') or $force or $empty or $regular){
    print "--> starting: $SEG_tool\n";
    seg();
}

if (($redo eq 'pfam') or $force or $empty or $regular){
    print "--> starting: $PFAM_tool\n";
    pfam();
}

if (($redo eq 'smart') or $force or $empty or $regular){
    print "--> starting: $SMART_tool\n";
    smart();
}

print "--> annotation finished.\n";

#### STOPPING TIME END ####
@time = localtime(time());
my $milsec2= time;
$milsec2 -= int($milsec2);
$milsec2 *= 100;
$milsec2 = int($milsec2);
if($time[1]<10){$time[1]="0".$time[1]}
$time[4]+=1;
if($time[4]<10){$time[4]="0".$time[4]}
my $END = $time[3]."/".$time[4]."/".($time[5]+1900)." ".$time[2].":".$time[1].":".$time[0].":".$milsec2;
print "tool start: $START\ntool end  : $END\n";
print "#####################\n\n";

########## SUBROUTINES ###########
sub extractAnnotation{
    my ($gene_id) = ($_[0]);
    # for all files
    my @annotationFiles = glob($annot."/*.xml");
    
    foreach(@annotationFiles){
        open (EXTRACT, "<".$_) or die ("ERROR: could not find file: $_. $!\n");
        my @content = <EXTRACT>;
        close EXTRACT;
        my ($base, $path, $suffix) = fileparse( $_, qr/\.[^.]*/ );
        mkpath($extract."/");
        open (WRITE_EX, ">".$extract."/".$base.$suffix) or die ("ERROR: could not create file: $extract/$qORp/$base$suffix. $!\n");

        for(my $ii=0;$ii<scalar(@content);$ii++){
            if ($content[$ii] =~ m/^<\?xml/){
                print WRITE_EX $content[$ii];
            }
            if ($content[$ii] =~ m/^<tool name=/){
                print WRITE_EX $content[$ii];
            }
            if($content[$ii] =~ m/<protein id=/){
                my @protein = split(/\"/,$content[$ii]);
                if($protein[1] eq $qORp){
                    my $jj = $ii;
                    while($content[$jj] !~ m/<\/protein>/){
                        print WRITE_EX $content[$jj];
                        $jj++;
                    }
                    if($content[$jj] =~ m/<\/protein>/){
                        print WRITE_EX $content[$jj];
                        last;
                    }
                }
            }
        }
        print WRITE_EX "</tool>";
        close WRITE_EX;
    }
}

=comment
 hashed %qids and %ids will be filled with gene ids from query and proteom fasta files
 TODO: adding versions in a smart way
 acquire seq length from input fasta
 fill hash with ids
=cut
sub seq_id_len {
    # creates hash lists of ids and seq length
    my $file = $_[0];
    my $id;
    
	open(INPUT, $file)   # open the input fasta file
	or die("ERROR: could not find or open file: $file. $!\n");
	my @content = <INPUT>;
	chomp(@content);
	for(my $i=0;$i<scalar(@content);$i++){  # loop through the content and find the IDs
		if($content[$i]=~/^>/){
			$content[$i] =~ s/\s+//g;               # substitute whitespaces
			$content[$i] =~ s/\n//g;                # substitute newlines
			my $CR = chr(13);			# define some weird character
			$content[$i] =~ s/\Q$CR//g;             # substitute some weird character
			$id = $content[$i];			# now $id is my clean identifier
			$id=~s/\>//;				# get rid of the ">" character
			$i++;					# go in the next line of your input file
			my $SQlen=0;				# set current length to zero
			while(defined($content[$i]) && !($content[$i]=~/^>/)){	# as long as the next sequence does not appear
				$content[$i] =~ s/ +//g;	# substitute whitespaces
				$SQlen+=length($content[$i]);	# get length of current line (stored in $content[$i]
				$i++;				# proceed with next line
			}$i--;
			$any_ids{$id}=$SQlen;			# set length in an hash: id -> length
		}
	}
	close INPUT;
}

sub smart{
    my $smartPATH    = $fasPath."/SMART";
    my $OutFilesPATH = $smartPATH."/output_files";
    my @content = ();
    chdir($smartPATH);

    require($SMART_tool);               # require: making subroutins from smart_scan.pl availible.
    my $outName = main_smart($fasta, $qORp);	# calling the main subroutine from smart_scan.pl due to previous call of require("smart_scan.pl")

    if(-e $OutFilesPATH."/".$outName){
        open(SMART, $OutFilesPATH."/".$outName);
        @content = <SMART>;
        close SMART;
    }else{
        print("ERROR: could not find or open $OutFilesPATH/smart.out. $!\n");
    }

    print "Deleting temporary output file...\n";
    my $delCommand = "rm -f ".$OutFilesPATH."/".$outName;
    system ($delCommand);

    chdir($fasPath);    

    open(OUT, ">".$dirOut."smart.xml") or print("ERROR: could not create output file smart.xml.\n");
    print OUT "<?xml version=\"1.0\"?>\n<tool name=\"SMART-DB\">\n";
    for(my $i=0;$i<@content;$i++){  # loop the content of smart_scan.out file and filter whats needed
        if($content[$i]=~/^>/){
            my %domainName;		# a hash that counts the instances of a clan in a sequence
            $content[$i]=~s/\>//;	# get rid of the ">" character ... done for all tools!!!
            $content[$i] =~ s/\s+//g;   # substitute whitespaces
            $content[$i] =~ s/\n//g;    # substitute newlines
            my $CR = chr(13);		# define some weird character
            $content[$i] =~ s/\Q$CR//g; # substitute some weird character
            
            my $length = check_id_length($content[$i]);
            print OUT "\t<protein id=\"" . $content[$i] . "\" length=\"" . $length . "\">\n";

            $i++;
            if($content[$i]=~/No hits detected that satisfy reporting thresholds/){
                #print $content[$i]."\n";
                print OUT "\t</protein>\n";
                next;
            }

            my %starts;
            my %ends;  # saves starts and ends to a clan / feature
            my @temp;
            my @temp_familyname;
            my %evalInfo;
            while(defined($content[$i]) && ($content[$i]=~/^\#/)){
                if(($content[$i]=~/family:/)){
                    #get family name first:
                    @temp_familyname = split(/family:/, $content[$i]);
                    @temp_familyname = split(/ +/, $temp_familyname[1]);  # $temp[1] is the family name now
                    # handle e-value family (feature type) bases
                    my @temp_fam_eval = split(/E-value:/, $content[$i]);
                    @temp_fam_eval = split(/ +/, $temp_fam_eval[1]);
                    chomp($temp_fam_eval[1]);
                    $evalInfo{$temp_familyname[1]} = $temp_fam_eval[1];

                    if(!defined($domainName{$temp_familyname[1]})){   # if clan (family?) name does not already exist in hash
                        $domainName{$temp_familyname[1]} = 0;
                    }
                }
                else{
                    while(defined($content[$i]) && ($content[$i]=~/^\#/) && !($content[$i]=~/family:/)){
                        my @position = split(/\|/, $content[$i]);   # $position[3] and $position[4] define start and end of domain, $position[11] defines the instance evalue
                        if(!defined($starts{$temp_familyname[1]})){$starts{$temp_familyname[1]} = $position[10]."@@".$position[3]}
                        else{$starts{$temp_familyname[1]} .= ", ".$position[10]."@@".$position[3]}
                        if(!defined($ends{$temp_familyname[1]})){$ends{$temp_familyname[1]} = "$position[4]"}
                        else{$ends{$temp_familyname[1]} .= ", $position[4]"}

                        my $InstCount = $domainName{$temp_familyname[1]};
                        $domainName{$temp_familyname[1]} = $InstCount+1;
                        $i++;
                    }$i--;
                }
                $i++;
            }
# take sorting into account
            my @keys = keys(%domainName);
            foreach(@keys){                
                print OUT "\t\t<feature type=\"" . $_ . "\" instance=\"" . $domainName{$_} . "\" clan=\"" . "---" . "\" evalue=\"" . $evalInfo{$_}. "\">\n";
                    my @temp_s = split(/, /, $starts{$_});
                    my @s;
                    my %e;
                    foreach(@temp_s){
                        my @e_s = split(/@@/,$_);
                        $e{$e_s[1]}=$e_s[0];
                        push(@s,$e_s[1]);
                    }
                    my @e = split(/, /, $ends{$_});
                    foreach(@s){
                        print OUT "\t\t\t<start inst_eval=\"$e{$_}\" start=\"$_\" end=\"".shift(@e)."\"/>\n";
                    }
                print OUT "\t\t</feature>\n";
            }
            print OUT "\t</protein>\n";
        }
    }
    print OUT "</tool>\n";
    close OUT;
}


sub pfam {
    my $PfamPATH = $fasPath."/Pfam";
    my $OutFilesPATH = $PfamPATH."/output_files";

    chdir($PfamPATH);

    require($PFAM_tool);
    my $outName = main($fasta,$qORp);	# call main from pfam_scan.pl (availible due to previous usage of require(pfam_scan.pl) -> pfam_scan_taxa.out

    open(PFAM, $OutFilesPATH."/".$outName) or print("ERROR: could not find or open $outName. $!\n $OutFilesPATH/$outName\n");
    my @content = <PFAM>;
    chomp(@content);
    close PFAM;

    print "Deleting temporary output file...\n";
    my $delCommand = "rm -f ".$OutFilesPATH."/".$outName;
    system($delCommand);

    chdir($fasPath);

    open(OUT, ">".$dirOut."pfam.xml") or print("ERROR: could not create output file pfam.xml.\n");
        
    print OUT "<?xml version=\"1.0\"?>\n<tool name=\"PfamDB\">\n";

    for(my $i=0;$i<@content;$i++){  # loop the content of pfam_scan.out file and filter whats needed
        if($content[$i]=~/^>/){
            my %pfamDomainName;		# a hash that counts the instances of a clan in a sequence
            $content[$i]=~s/\>//;	# get rid of the ">" character ... done for all tools!!!
            $content[$i] =~ s/\s+//g;   # substitute whitespaces
            $content[$i] =~ s/\n//g;    # substitute newlines
            my $CR = chr(13);		# define some weird character
            $content[$i] =~ s/\Q$CR//g; # substitute some weird character
            
            my $length = check_id_length($content[$i]);
            print OUT "\t<protein id=\"" . $content[$i] . "\" length=\"" . $length . "\">\n";

            $i++;
            if($content[$i]=~/No hits detected that satisfy reporting thresholds/){
                #print $content[$i]."\n";
                print OUT "\t</protein>\n";
                next;
            }

            my %starts;
            my %ends;  # saves starts and ends to a clan / feature
            my @temp;
            my @temp_familyname;
            my %clanInfo;
            my %evalInfo;
            while(defined($content[$i]) && ($content[$i]=~/^\#/)){
                if(($content[$i]=~/family:/)){
                    #get family name first:
                    @temp_familyname = split(/family:/, $content[$i]);
                    @temp_familyname = split(/ +/, $temp_familyname[1]);  # $temp[1] is the family name now
                    # handle e-value family (feature type) bases
                    my @temp_fam_eval = split(/E-value:/, $content[$i]);
                    @temp_fam_eval = split(/ +/, $temp_fam_eval[1]);
                    $evalInfo{$temp_familyname[1]} = $temp_fam_eval[1];
                    #handling clan info
                    if(!($content[$i]=~/clan: ---/)){		#clan info is av.
                        @temp = split(/clan:/, $content[$i]);
                        @temp = split(/ +/, $temp[1]);  # $temp[1] is the clan name now

                        if(!defined($clanInfo{$temp_familyname[1]})){ #clan detected the first time
                            $clanInfo{$temp_familyname[1]} = $temp[1];
                        }
			
                    }elsif($content[$i]=~/clan: ---/){ #no clan info av. (family name will be used as "clan" to 
                        #set family name as clan information
                        if(!defined($clanInfo{$temp_familyname[1]})){ #clan detected the first time
                            $clanInfo{$temp_familyname[1]} = $temp_familyname[1];
                        }
			        }

                    if(!defined($pfamDomainName{$temp_familyname[1]})){   # if clan (family?) name does not already exist in hash
                        $pfamDomainName{$temp_familyname[1]} = 0;
                    }
                }
                else{
                    while(defined($content[$i]) && ($content[$i]=~/^\#/) && !($content[$i]=~/family:/)){
                        #print "$i: ".$content[$i]."\n";
                        my @position = split(/\|/, $content[$i]);   # $position[3] and $position[4] define start and end of domain, $position[11] defines the instance evalue
                        if(!defined($starts{$temp_familyname[1]})){$starts{$temp_familyname[1]} = $position[10]."@@".$position[3]}
                        else{$starts{$temp_familyname[1]} .= ", ".$position[10]."@@".$position[3]}
                        if(!defined($ends{$temp_familyname[1]})){$ends{$temp_familyname[1]} = "$position[4]"}
                        else{$ends{$temp_familyname[1]} .= ", $position[4]"}

                        my $InstCount = $pfamDomainName{$temp_familyname[1]};
                        $pfamDomainName{$temp_familyname[1]} = $InstCount+1;
                        $i++;
                    }$i--;
                }
                $i++;
            }
# take sorting into account
            my @keys = keys(%pfamDomainName);
            foreach(@keys){                
                #print OUT "\t\t<feature type=\"$_\" instance=\"$pfamDomainName{$_}\">\n";
                print OUT "\t\t<feature type=\"" . $_ . "\" instance=\"" . $pfamDomainName{$_} . "\" clan=\"" . $clanInfo{$_} . "\" evalue=\"" . $evalInfo{$_}. "\">\n";
                    my @temp_s = split(/, /, $starts{$_});
                    my @s;
                    my %e;
                    foreach(@temp_s){
                        my @e_s = split(/@@/,$_);
                        $e{$e_s[1]}=$e_s[0];
                        push(@s,$e_s[1]);
                    }
                    my @e = split(/, /, $ends{$_});
                    foreach(@s){
                        print OUT "\t\t\t<start inst_eval=\"$e{$_}\" start=\"$_\" end=\"".shift(@e)."\"/>\n";
                        #print OUT "\t\t\t<end end=\"".shift(@e)."\"/>\n";
                    }
                print OUT "\t\t</feature>\n";
            }
            print OUT "\t</protein>\n";
        }
    }
    print OUT "</tool>\n";
    close OUT;
}



sub tempfile {
# function returns an array with ID and sequence of each protein in the input file
    my $file = $fasta;
    my @return;
    
    open(FILE, $file)
        or print("ERROR: could not find input file.\n");
    my @content = <FILE>;
    chomp(@content);
    close FILE;

    for(my $i=0;$i<@content;$i++){  # loop through input file
        if(defined($content[$i]) && $content[$i]=~/^\>/){   # if an id was found
            push(@return, $content[$i]);    # save the id
            $i++;                           # go to next line
            my $tempSeq="";                    # temp save for the sequence
            while(defined($content[$i]) && !($content[$i]=~/^\>/)){ # while this line is not a new ID
                $tempSeq.=$content[$i];    # add seq to temp save
                $i++;
            }$i--;
            push(@return, $tempSeq);        # save the sequence            
        }
    }

    return(@return);    # @return (ID, seq, ID2, seq2...)
}

sub seg {
# this function calls seg (compiled c++) which identifies low complexity regions
# on a protein sequence;
    my $segPath = $fasPath."/SEG";

    my $command = "\"$segPath/$SEG_tool\" \"$fasta\" -l -n -p |";

    open(RESULT, $command)
        or print("ERROR: could not start seg-tool.\n");

    my @result = <RESULT>;

    chomp(@result);
    close RESULT;

    open(OUT, ">".$dirOut."seg.xml")  # create output file, overwrite existing one
        or print("ERROR: could not write output file for seg.\n");

    print OUT "<?xml version=\"1.0\"?>\n<tool name=\"seg\">\n";

    for(my $i=0;$i<@result;$i++){
        if($result[$i]=~/^>/){   # find ID
            my @output;
            if($i>0){print OUT "\n"}
            $result[$i]=~s/\>//;
            $result[$i] =~ s/\s+//g;        # substitute whitespaces
            $result[$i] =~ s/\n//g;         # substitute newlines
            my $CR = chr(13);               # define some weird character
            $result[$i] =~ s/\Q$CR//g;      # substitute some weird character

            my $length = check_id_length($result[$i]);
            print OUT "\t<protein id=\"" . $result[$i] . "\" length=\"" . $length . "\">\n";            
            
            $i++;

            my $lcrCounter = 0;
            while(defined($result[$i]) && !($result[$i]=~/^>/)){  # loop through the results to this ID
                if(!($result[$i] eq "")){ # if line is not empty
                    $result[$i]=~s/^ +//;   # delete whitespaces at line beginning
                    my @temp = split(/ +/, $result[$i]);
                    if(defined($temp[0]) && $temp[0]=~/[a-z]/ && defined($temp[1]) && $temp[1]=~/[0-9]/){
                        $lcrCounter++;
                        my @tempPos = split(/-/, $temp[1]);  
                        push(@output, "\n\t\t\t<start start=\"$tempPos[0]\"/>\n\t\t\t<end end=\"$tempPos[1]\"/>");
                    }
                }
                $i++;
            }$i--;
            print OUT "\t\t<feature type=\"low complexity regions\" instance=\"$lcrCounter\">";
            foreach(@output){print OUT $_;}  # print into output file
            print OUT "\n\t\t</feature>";
            print OUT "\n\t</protein>";
        }
    }print OUT "\n</tool>";

    close OUT;
}

sub signalp {
# this function calls the signalp.pl script to find signal peptides in protein 
# output will be parsed into xml format
    my $SigPpath = $fasPath."/SignalP";
    my @idlist;

    open(INPUT, "$fasta")
        or print("ERROR: function: signalp() / opening input file. $!\n");
    my @temp = <INPUT>;
    foreach(@temp){
        if($_=~/\>/){$_=~s/\>//;push(@idlist, $_)}        
    }
    chomp(@idlist);
    close INPUT;

    chdir($SigPpath);
    my @result	= `perl $SigPpath/$SIGNALP_tool -t euk \"$fasta\"`;
	chomp(@result);

    chdir($fasPath);

    open(OUT, ">".$dirOut."signalp.xml")  # create output file, overwrite existing one
        or print("ERROR: could not write output file for signalp. $!\n");

    print OUT "<?xml version=\"1.0\"?>\n<tool name=\"signalp\">\n";

    for(my $i=0;$i<scalar(@result);$i++){
        if(!($result[$i]=~/^#/)){
            my @line = split(/\s+/, $result[$i]);
            foreach(@idlist){
                my $tempID = $_;
                if ($tempID eq $line[0]){
                    $tempID=~s/\>//;
                    $tempID =~ s/\s+//g;         # substitute whitespaces
                    $tempID =~ s/\n//g;         # substitute newlines
                    my $CR = chr(13);			# define some weird character
                    $tempID =~ s/\Q$CR//g;      # substitute some weird character

                    my $length = check_id_length($tempID);
                    if ($line[9] eq "Y"){
                        #end of feature = Ymaxpos -1 (as is signalp.pl version 4.1)
                        my $end = $line[4]-1;
                        #yes case of SigP
                        print OUT "\t<protein id=\"" . $tempID . "\" length=\"" . $length . "\">\n";
                        print OUT "\t\t<feature type=\"SIGNAL\" instance=\"1\">\n";
                        print OUT "\t\t\t<start start=\"1\"/>\n\t\t\t<end end=\"$end\"/>\n";
                        print OUT "\t\t</feature>\n";
                        print OUT "\t</protein>\n";

                    }elsif($line[9] eq "N"){
                        #no case of SigP
                        print OUT "\t<protein id=\"" . $tempID . "\" length=\"" . $length . "\">\n";
                        print OUT "\t\t<feature type=\"SIGNAL\" instance=\"0\">\n";
                        print OUT "\t\t</feature>\n";
                        print OUT "\t</protein>\n";

                    }else{
                            #signalP output might be strange
                    }
                }   
            }
        }
    }
    print OUT "</tool>";
    close OUT;
}

sub coils {
# this function calls coils-tool which finds coiled coil regions in the protein sequence
# and then creates an output with the results. legende: c - coil number, s - start point (aa)
# of this coil, e - end point (aa) of this coil, l - length of this coil, tc - total count
# of coils in this protein sequence
    my $COILSpath=$fasPath."/COILS2";
    chdir($COILSpath);  # COILS2-tool has to be called from its directory, otheriwse it cant find its binaries

    my $command = "\"$COILSpath/$COILS_tool\" -f -win 21 < \"$fasta\" |";

    open(RESULT, $command) # start coils2 and get result in filehandler
        or print ("ERROR: could not start coils.\n");
    my @result = <RESULT>;
    chomp(@result);
    close RESULT;
    
    chdir($fasPath);
    
    open(OUT, ">".$dirOut."coils.xml")  # create output file, overwrite existing one
        or print("ERROR: could not write output file for coils.\n");

    print OUT "<?xml version=\"1.0\"?>\n<tool name=\"COILS2\">\n";

    for(my $i=0;$i<@result;$i++){   # now writing in final output file: "coils.out"
        if($i>0){$i--}
        if($result[$i]=~/\>/){  # find ID
            $result[$i]=~s/\>//;
            $result[$i]=~s/\s+//g;
            $result[$i] =~ s/\n//g;     # substitute newlines
            my $CR = chr(13);			# define some weird character
            $result[$i] =~ s/\Q$CR//g;  # substitute some weird character

            my $length = check_id_length($result[$i]);
            print OUT "\t<protein id=\"" . $result[$i] . "\" length=\"" . $length . "\">\n";
            
            $i++;   # go to next line where the sequence should start

            my $sequence;
            while(defined($result[$i]) && !($result[$i]=~/\>/)){ # making a one-line-sequence
                $sequence.=$result[$i];
                $i++;
            }

            my @temp = split(/x/, $sequence);   # split output at cc positions
            my $coilCounter=0;  # counts the number of coils
            my @pos;    # save all starts and ends to print
            if(defined($temp[1])){  # do the rest only if at least one coiled coil position is found
                my $position=0; # temp save for start and end positions of coils
                for(my $i=0;$i<@temp;$i++){
                    if($i>0){$i-=1}
                    if(!($temp[$i] eq "") or (($i==0) && ($temp[$i] eq ""))){ # if line is not empty or is empty but its the first line
                        $coilCounter++;
                        $position+=length($temp[$i])+1;   # start position of coil is length of this line (plus end position from before)
                        push(@pos, $position);
                        $i++;   # go to next line

                        if(defined($temp[$i]) && $temp[$i] eq ""){    # now this line should be empty
                            my $coilLEN=1;      # length of coil counter
                            while($temp[$i] eq ""){ # as long as there are empty lines
                                $coilLEN++; # increase the coil length counter
                                $i++;   # and go to next line
                            }
                            $position+=($coilLEN-1);
                            push(@pos, $position);
                        }
                    }
                }
            }
### OUTPUT :
            print OUT "\t\t<feature type=\"coiled_coil\" instance=\"$coilCounter\">\n";
            my $x=0;
            foreach(@pos){
                if($x==0){print OUT "\t\t\t<start start=\"$_\"/>\n";$x++;}
                elsif($x==1){print OUT "\t\t\t<end end=\"$_\"/>\n";$x--;}
            }
            print OUT "\t\t</feature>\n";
        }print OUT "\t</protein>\n";
    }
    print OUT "</tool>";
    close OUTPUT;
}

sub TMHMM20 {   
# create TMHMM2.0 output - finds Transmembrane regions, doesnt say anything 
# about proteins without TR-regions
    my $TMHMMpath=$fasPath."/TMHMM"; # path to compiled TMHMM file (decodeanhmm)

    my $command = "cat \"$fasta\" | \"$TMHMMpath/$TMHMM_tool\" -f \"$TMHMMpath/TMHMM2.0.options\" -modelfile \"$TMHMMpath/TMHMM2.0.model\" |";

    open(RESULT, $command) # start decodeanhmm and get result in filehandler
        or print ("ERROR: could not start TMHMM2.0(decodeanhmm).\n");
    my @result = <RESULT>;
    chomp(@result);
    close RESULT;

    open(OUT, ">".$dirOut."tmhmm.xml")  # create output file, overwrite existing one
        or print("ERROR: could not write output file for tmhmm.\n");

    print OUT "<?xml version=\"1.0\"?>\n<tool name=\"TMHMM2.0c\">\n";

    for(my $i=0;$i<@result;$i++) {  # go through the results, filter out whats needed and create output file
        my $TMcount=0;  # counts the number of TR regions
        if($result[$i]=~/\>/){  # find ID
            $result[$i] =~ s/\>//;
            $result[$i] =~ s/\s+//g;        # substitute whitespaces
            $result[$i] =~ s/\n//g;         # substitute newlines
            my $CR = chr(13);               # define some weird character
            $result[$i] =~ s/\Q$CR//g;      # substitute some weird character

            my $length = check_id_length($result[$i]);
            print OUT "\t<protein id=\"" . $result[$i] . "\" length=\"" . $length . "\">\n"; 
            
            $i++;
            if($result[$i]=~/\%/){  # find line with information of interest
                my @temp = split(/\:/, $result[$i]);
                @temp = split(/\, /, $temp[1]);
                $temp[0]=~s/ //;
                chomp(@temp);
                foreach(@temp){
                    if($_=~/M/){$TMcount++}
                    $_=~s/^O/o/;
                    #$_.="\n";
                }
### OUTPUT :
                print OUT "\t\t<feature type=\"transmembrane\" instance=\"$TMcount\">\n";
                foreach(@temp){
                    my @pos;
                    if($_=~/M/){
                        @pos = split(/ +/, $_);
                        print OUT "\t\t\t<start start=\"$pos[1]\"/>\n";
                        print OUT "\t\t\t<end end=\"$pos[2]\"/>\n";
                    }
                }
                print OUT "\t\t</feature>\n";
            }print OUT "\t</protein>\n";
        }
    }
    print OUT "</tool>\n";
    close OUTPUT;
}
#check_id_length:   acquire length from precalculated hash
#check for match/agreement/concurrence of identifier
sub check_id_length{
	my $current_id = $_[0];
	
	if (defined($any_ids{$current_id})){
            return $any_ids{$current_id};
	}else{
            die ( "\nError: " . $current_id . " is a strange identifier ... exiting. Annotation went wrong. Please check your sequence IDs/header.\n\n");
	}
}

#TODO: infile documentation 
# create CAST ouput - finds regions enriched for a paritcular AA, threshold is important!
sub CASTing {   
	my $threshold=50; #default
	my $castpath=$fasPath."/CAST"; # path to compiled cast file 

	unless(-d "$castpath/tmp"){
		system("mkdir $castpath/tmp");
	}

	my $result; 
	
	# open protein fasta input
	open(READ,$fasta) || die "Cannot open $fasta!\n";
	my @fas = <READ>;
	close (READ);

	# get speciesName
	my @fastaTMP = split(/\//,$fasta);
	my @speciesNameTMP = split(/\./,$fastaTMP[@fastaTMP-1]);
	my $specName = $speciesNameTMP[0];

	# split into multiple sequences
	my $fas = join("",@fas);
	my @allSeq = split(">",$fas);

	# run CAST for each sequence and save output into @result
        my $pid = $$;
        my $tempfasta = $castpath."/tmp/".$specName."_".$pid."tmp.fa";
	foreach my $seq (@allSeq){
		if(length($seq)>2){
			chomp($seq);
			open(TMP,">$tempfasta");
			print TMP ">",$seq;
			close TMP;
		
			my $castOut =  qx($castpath/$CAST_tool $tempfasta -verbose -thr $threshold);
			chomp($castOut);
			$result .= $castOut."\n";
		}
	}

	chomp($result);
	my @result = split(/\n/,$result);

	# write output into XML format
	open(OUT, ">".$dirOut."cast.xml")   # create output file, overwriting if one exists
	    or print ("ERROR: could not write output file for CAST. $!\n");

	print OUT "<?xml version=\"1.0\"?>\n<tool name=\"CAST\">\n";
        ## results will be written in xml output file
	for(my $i=0;$i<@result;$i++) {
	    my @value;
	    if($result[$i]=~/\>/){  
                chomp($result[$i]);
                $result[$i]=~s/\>//;
                $result[$i] =~ s/\s+//g;	# substitute whitespaces
                my $CR = chr(13);		# define some weird character
                $result[$i] =~ s/\Q$CR//g;	# substitute some weird character

                my $length = check_id_length($result[$i]);
                print OUT "\t<protein id=\"" . $result[$i] . "\" length=\"" . $length . "\">\n";   # print id
                $i++;
                while($result[$i] =~ m/(.*) region from (\d*) to (\d*) corrected with score (\d*)/){
                    @value = split(/ +/, $result[$i]);
                    chomp(@value);
                    print OUT "\t\t<feature type=\"$value[0]\" instance=\"1\">\n";
                    print OUT "\t\t\t<start start=\"$value[3]\"/>\n";
                    print OUT "\t\t\t<end end=\"$value[5]\"/>\n";
                    print OUT "\t\t</feature>\n";
                    $i++;
                }
                print OUT "\t</protein>\n";
	    }
	}
    print OUT "</tool>\n";
	close OUTPUT;
	close RESULT;
	system("rm -f $tempfasta");
}
################
sub getHelpMessage{
    my $message=
"You are running $0 version $version\n\n$0 annotates a given sequence file. Please define your sequence file (FASTA format),
the path for the output directory and the name for your annotation output.\n

Examples:\n
Run $0 on a fasta file (multifasta possible):
$0 -fasta=/path/to/your/fastafile -path=/path/to/directory -name=givenname\n
Extract from existing annotations:
$0 -path=/path/to/your/existing/annotationdir -name=genename -extract=/path/to/directory/

Options:
 -h
            print this message,
 -fasta=<>
            specify the sequence file,
 -path=<>
            specify the location where your annotation files will be stored,
 -name=<>
            specify a name for query.

ADDITIONAL OPTIONS

 -force
            set this flag if you want to force annotations (override),
 -extract=<>
            specify a path to the location where you want the extracted annotations to be stored.
 -redo=<>
            specifiy for which feature database you want to re-annotate the sequence file. [cast|coils|seg|pfam|signalp|smart|tmhmm] (Only one selection possible)\n";      
            
return $message;
}
