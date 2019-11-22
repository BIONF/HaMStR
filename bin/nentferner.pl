#!/home/vinh/anaconda3/envs/hamstr/bin/perl
use lib '/home/vinh/programs/HaMStR/lib';
use Filehandler;
use strict;
use Getopt::Long;

# PROGRAM NAME: NENTFERNER.PL

# AUTHOR: INGO EBERSBERGER, ebersber@eva.mpg.de
# DATE: 
# DESCRIPTION: THE PROGRAM ALL REMOVES NEWLINES FROM THE INPUT TEXT

# DATE LAST MODIFIED: 22/01/2001; 12.02.2004
## Last modified: 09.01.2014
## added the option to define an outfile including an outpath
#################### START MAIN PROGRAM #####################


if (-e "nentferner.out") {
    `rm -f nentferner.out`;
}
my $pid = $$;
my $p;
## parse the command line
my $infile;
my $space = '';
my $sep = '';
my $path2infile;
my $path2outfile;
my $ls;
my $rfq;
my $outname;
my $help;
my $ln;
my $replace;
my $offset = 1;
GetOptions ("h" => \$help,
	    "in=s" => \$infile, 
	    "out=s"=> \$outname,
	    "space=s" => \$space,
	    "sep=s" => \$sep,
	    "leading_space" => \$ls,
	    "reformat_qual" => \$rfq,
	    "line_numbers" => \$ln,
	    "offset=s" => \$offset,
	    "replace" => \$replace);
 
###############
## help
if (defined $help) {
    die "options:
-in=<path2infile/infile>\n
-out=<path2outfile/outfile>\n
-space=<yes||no>\tdefault:no\n
-sep=<tab||newline||no>\tdefault:newline. Set value to 'no' if you just want to join the lines\n
-leading_space: set this flag to remove one blank from the beginning of each line in the file\n
-reformat_qual: set this flag to reformat phred qual files such that every qual value occupies two bytes\n
-line_numbers\tChoose this option if you want to have the line number added to the output. Works best with -sep=tab. 
-offset=<>\tGive the number to start the line numbering with. Default is 1
-replace\tSet this flag if you want to remove the newlines in the original file\n";
}

my $seq;
my $crunch = 0;
my @head;
my $intcount;
my $join;
## interpretation of command line values
if ($infile =~ /\//) {
    ($path2infile, $infile) = $infile =~ /(.*)\/(.*)?$/;
}
if ($outname =~ /\//) {
    ($path2outfile, $outname) = $outname =~ /(.*)\/(.*)?$/;
}

if (!(defined $path2infile)) {
    $path2infile = '.';
}
if (!(defined $path2outfile)) {
    $path2outfile = $path2infile;
}

#
if ($space eq 'yes') {
	$space = ' ';
}
else {
	$space = '';
	print "\nnewlines will be removed without replacement\n";
}
#
if ($ln and !$sep) {
	$sep = 'tab';
}
if ($sep eq 'tab') {
    $join = "\t";
    $sep = "\t";
}
elsif ($sep eq 'no') {
    $join = '';
    $sep = '';
}
else {
    $join = ' ';
    $sep = "\n";
}
#i####### if reformat_qual is set, automatically set $sep to ' ' and $ls to 1
if ($rfq) {
	$space = ' ';
	$ls = 1;
}
if (!(defined $outname)) {
    $outname = 'nentferner.out';
}

die "The file $path2infile/$infile does not exist!\n" unless -e "$path2infile/$infile";
tie (*IN, Filehandler::, "$infile", "$path2infile", "\n");
tie (*OUT, Filehandler::, "$outname", ">$path2outfile", "\n");
while (readline (IN)) {
    
    if ($_ =~ />/) {
	## the fasta header
	if ($crunch == 1) {
	    if (defined $rfq) {
            $seq =~ s/(\s{1}\d{1}\s{1})/ $1/g;
            $seq =~ s/(\s{1}\d{1})(\s{1}\d{1}\s{1})/$1 $2/g;
	    $seq =~ s/^(\d{1}\s{1})/ $1/;
	    $seq =~ s/(\s{1}\d{1})$/ $1/;
	    }
	    my $outline = (join "$join", @head) . "$sep" . $seq . "\n";
	    print OUT $outline;
	}
	chomp $_;
	my $head = $_;
#	$head =~ s/>.*\|//;
	@head = split /\s{1,}/, $head;
	if ($ln) {
		$head[0] =~ s/>//;
		my $no=$intcount + $offset;
		@head=($no,@head);
	}
	$seq = '';
	$crunch = 1;
	$intcount ++;
	if ($intcount%1000 == 0) {
	    print "$intcount lines processed\n";
	}
    }
    else {
	if (defined $ls) {
		$_ =~ s/^\s{1}//;
	}
	$_ =~ s/\s{1,}$/$space/g;
	$seq .= $_;
    }
}
if (defined $rfq) {
            $seq =~ s/(\s{1}\d{1}\s{1})/ $1/g;
            $seq =~ s/(\s{1}\d{1})(\s{1}\d{1}\s{1})/$1 $2/g;
            $seq =~ s/^(\d{1}\s{1})/ $1/;
            $seq =~ s/(\s{1}\d{1})$/ $1/;
            }

my $outline = (join "$join", @head) . "$sep" . $seq . "\n";
print OUT $outline;
close (OUT);
close (IN);
if ($replace) {
	`mv $path2outfile/$outname $path2infile/$infile`; 
}
exit;
