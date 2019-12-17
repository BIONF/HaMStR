#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

flag=0
root=1
### check grep, sed and wget availability
echo "-------------------------------------"
echo "Checking .bash_profile/.bashrc, grep, sed/gsed and wget availability..."
grepprog='grep'
sedprog='sed'
wgetprog='wget'
bashFile='.bashrc'
if [ "$sys" == "Darwin" ]; then
    if [ -z "$(which brew)" ]; then
        echo "Please install homebrew to install dependencies tools and libraries!"
        echo "Check https://brew.sh"
        exit
    fi
    sedprog='gsed'
	grepprog='ggrep'
	shell=$(echo $SHELL)
	if [ $shell == "/bin/zsh" ]; then
    	bashFile='.zshrc'
	else
		bashFile='.bash_profile'
	fi
else
    if [ "$EUID" -ne 0 ]; then
        echo "You are not running this setup as root."
        read -p "Press enter to continue, but some missing tools/libraries will not be installed!"
        root=0
    fi
fi

if [ -z "$(which $sedprog)" ]; then
    if [ "$sys" == "Darwin" ]; then
        brew install gnu-sed
    fi
fi

if [ -z "$(which $grepprog)" ]; then
    if [ "$sys" == "Darwin" ]; then
        brew install grep
    fi
fi

if [ -z "$(which $wgetprog)" ]; then
    if [ "$sys" == "Darwin" ]; then
        brew install wget
    fi
fi

if ! [ -f ~/$bashFile ]; then
    touch ~/$bashFile
fi
if [ "$flag" == 1 ]; then exit 1; fi
echo "done!"

### check dependencies
echo "-------------------------------------"
echo "Installing dependencies..."

dependenciesUbuntu=(
  build-essential # for make
  curl
  r-base # for Statistics::R
  wise
  hmmer # hmmer (for both hmmsearch and hmmbuild)
  clustalw
  mafft
  muscle
  blast2 # blast
  ncbi-blast+
  libdbi-perl
  libipc-run-perl
  perl-doc
  locales
  lib32ncurses5
  lib32z1
)

dependenciesMac=(
  brewsci/bio/genewise
  hmmer # hmmer (for both hmmsearch and hmmbuild)
  brewsci/bio/clustal-w
  mafft
  brewsci/bio/muscle
  blast
)

if [ "$sys" == "Darwin" ]; then
  for i in "${dependenciesMac[@]}"; do
  	echo $i
  	brew install $i
  done
  if [ -z "$(grep clustalw ~/$bashFile)" ]; then
      echo "alias clustalw='clustalw2'" >> ~/$bashFile
  fi
else
    if [ $root == 1 ]; then
        sudo apt-get update -y
        for i in "${dependenciesUbuntu[@]}"; do
        	echo $i
        	sudo apt-get install -y -qq $i > /dev/null
        done
    fi
fi

dependencies=(
  genewise
  hmmsearch
  hmmbuild
  mafft
  muscle
  blastn
)

for i in "${dependencies[@]}"; do
  if [ -z "$(which $i)" ]; then
    echo "$i not found / cannot be automatically installed. Please install it and run this setup again!"
    flag=1
  fi
done
if [ "$flag" == 1 ]; then exit 1; fi

wisePath=$(which "genewise")
if [ -z "$(grep WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
    echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
fi

echo "Installing Perl modules..."
perlModules=(
  DBI
  DB_File
  File::Copy
  File::Path
  File::Basename
  File::Which
  List::Util
  Parallel::ForkManager
  POSIX
  XML::SAX
  XML::NamespaceSupport
  XML::Parser
  Getopt::Long
  IO::Handle
  IPC::Run
  Statistics::R
  Term::Cap
  Time::HiRes
  Bio::AlignIO
  Bio::Align::ProteinStatistics
  Bio::DB::Taxonomy
  Bio::SearchIO
  Bio::SearchIO::blastxml
  Bio::Search::Hit::BlastHit
  Bio::Seq
  Bio::SeqIO
  Bio::SeqUtils
  Bio::Tree::Tree
  Bio::Tools::Run::StandAloneBlast
)

if [ "$sys" == "Darwin" ]; then
    if [ -z "$(which cpanm)" ]; then
      curl -L http://cpanmin.us | perl - --sudo App::cpanminus
    fi

    for i in "${perlModules[@]}"; do
      msg=$((perldoc -l $i) 2>&1)
      if [[ "$(echo $msg)" == *"No documentation"* ]]; then
        sudo cpanm ${i} --quiet --force
      fi
    done
else
    if [ $root == 1 ]; then
        if [ -z "$(which cpanm)" ]; then
          curl -L http://cpanmin.us | perl - --sudo App::cpanminus
        fi

        for i in "${perlModules[@]}"; do
          msg=$((perldoc -l $i) 2>&1)
          if [[ "$(echo $msg)" == *"No documentation"* ]]; then
            sudo cpanm ${i} --quiet --force
          fi
        done
    fi
fi
echo "done!"

### prepare folders
echo "-------------------------------------"
echo "Preparing folders..."
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR/..
CURRENT=$(pwd)

# create required folders
folders=(
  blast_dir
  core_orthologs
  genome_dir
  weight_dir
  taxonomy
  output
  tmp
  "bin/aligner"
)

for i in "${folders[@]}"; do
  echo "$i"
  if [ ! -d $i ]; then mkdir $i; fi
done
echo "done!"

### download tools
echo "-------------------------------------"
echo "Downloading and installing annotation tools/databases:"

fasta36="yes"
if [ -z "$(which fasta36)" ]; then
  fasta36="no"
  fasta36v="fasta-36.3.8h"
  if ! [ -f "bin/aligner/bin/fasta36" ]; then
	  echo "fasta-36"
	  wget "http://faculty.virginia.edu/wrpearson/fasta/fasta36/${fasta36v}.tar.gz"
	  tar xfv $fasta36v.tar.gz
	  rm "${fasta36v}.tar.gz"
	  mv $fasta36v/* bin/aligner/
	  rm -rf $fasta36v
	  cd "bin/aligner/src"
	  if [ $sys=="Linux" ]; then
	    make -f ../make/Makefile.linux64_sse2 all
	  elif [ $sys=="Darwin" ]; then
	    make -f ../make/Makefile.os_x86_64 all
	  fi
  fi
  if [ -z "$(grep PATH=$CURRENT/bin/aligner/bin ~/$bashFile)" ]; then
	  echo "export PATH=$CURRENT/bin/aligner/bin:\$PATH" >> ~/$bashFile
  fi
fi
cd $CURRENT
if [ -z "$(which fasta36)" ]; then
	if ! [ -f "$CURRENT/bin/aligner/bin/fasta36" ]; then
		echo "fasta36 tool could not be found in $CURRENT/bin/aligner/. Please check again!"
		exit
	fi
fi

cd "taxonomy"
if ! [ -f "nodes" ]; then
    wget "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tar xfv taxdump.tar.gz
    rm taxdump.tar.gz
    echo "Taxonomy database indexing. It can take a while, please wait..."
    perl $CURRENT/bin/indexTaxonomy.pl $CURRENT/taxonomy
    rm *.dmp
    rm gc.prt
    rm readme.txt
fi
cd $CURRENT
if ! [ -f "$CURRENT/taxonomy/nodes" ]; then
	echo "Error while indexing NCBI taxonomy database! Please check $CURRENT/taxonomy/ folder and run this setup again!"
	exit
fi

cd "bin"
if [ -z "$(which greedyFAS)" ]; then
# 	# if ! [ -f "fas/Pfam/Pfam-hmms/Pfam-A.hmm"]; then
	    echo "FAS"
	    wget https://github.com/BIONF/FAS/archive/master.tar.gz
	    tar xfv master.tar.gz
# 		rm -rf fas
	    mv FAS-master fas
	    rm master.tar.gz
        pip install $CURRENT/bin/fas
        if [ -z "$(which annoFAS)" ]; then
            echo "Installation of FAS failed! Please try again!"
            exit
        else
            annoFAS --fasta test.fa --path $CURRENT --name q --prepare 1 --annoPath $CURRENT/bin/fas
        fi
# 	    # chmod 755 fas/config/setup.sh
# 	    # fas/config/setup.sh
# 	# fi
fi
cd $CURRENT
if ! [ -f "$CURRENT/bin/fas/Pfam/Pfam-hmms/Pfam-A.hmm" ]; then
	echo "Annotation tools are missing! Please install them again using"
    echo "annoFAS --fasta pseudo.fa --path pseudo.path --name q --prepare 1 --annoPath $CURRENT/fas"
    exit
fi
echo "done!"

### download data
echo "-------------------------------------"
echo "Getting pre-calculated data"

if ! [ "$(ls -A $CURRENT/genome_dir)" ]; then
	echo "Processing $CURRENT ..."
	if [ ! -f $CURRENT/data_HaMStR.tar ]; then
		echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
		wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar
	else
		CHECKSUM=$(cksum data_HaMStR.tar)
		echo "Checksum: $CHECKSUM"
		if ! [ "$CHECKSUM" == "4100986910 5840435200 data_HaMStR.tar" ]; then
    		  rm $CURRENT/data_HaMStR.tar
    		  echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
      		  wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar
    	fi
    fi

	if [ ! -f $CURRENT/data_HaMStR.tar ]; then
	  echo "File data_HaMStR.tar not found! Please try to download again from"
	  echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
	  exit
	fi

	CHECKSUM=$(cksum data_HaMStR.tar)
	if [ "$CHECKSUM" == "4100986910 5840435200 data_HaMStR.tar" ]; then
	  echo "Extracting archive data_HaMStR.tar"
	  tar xfv $CURRENT/data_HaMStR.tar
	  rm $CURRENT/data_HaMStR.tar
	  echo "Archive data_HaMStR.tar extracted into $CURRENT"
	  if [ ! -d $CURRENT/data_HaMStR ]; then
		echo "Directory $CURRENT/data_HaMStR not found!"
	  else
		printf "\nMoving gene sets ...\n--------------------\n"
		rsync -rva data_HaMStR/genome_dir/* $CURRENT/genome_dir
		printf "\nMoving blast databases ...\n--------------------------\n"
		rsync -rva data_HaMStR/blast_dir/* $CURRENT/blast_dir
		printf "\nMoving annotations ...\n----------------------\n"
		rsync -rva data_HaMStR/weight_dir/* $CURRENT/weight_dir
		printf "\nRemoving duplicated data. Please wait.\n------------------------------------\n"
		rm -rf $CURRENT/data_HaMStR
		printf "\nDone! Data should be in place to run HaMStR.\n"
	  fi
	else
	  echo "Something went wrong with the download. Checksum does not match."
	  echo "Please try to download again from"
	  echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
	  echo "Please put it into $CURRENT folder and run this setup again!"
	  exit
	fi
fi

### add paths to bash profile file
echo "-------------------------------------"
echo "Adding paths to ~/$bashFile"

if [ -z "$(grep PATH=$CURRENT/bin:\$PATH ~/$bashFile)" ]; then
	echo "export PATH=$CURRENT/bin:\$PATH" >> ~/$bashFile
fi

wisePath=$(which "genewise")
if [ -z "$(grep WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
    echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
fi
echo "done!"

### adapt paths in hamstr scripts
echo "-------------------------------------"
echo "Adapting paths in hamstr scripts"
# update the sed and grep commands
$sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/hamstr.pl;
$sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/hamstr.pl;
$sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/oneSeq.pl;
$sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/oneSeq.pl;

# localize the perl installation
path2perl=`which perl`
echo "path to perl: $path2perl"
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/hamstr.pl;
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/nentferner.pl;
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/translate.pl;
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/oneSeq.pl;

# get lib path
path2dir=$CURRENT
echo "path to lib: $path2dir/lib"
$sedprog -i -e "s|use lib.*lib\(.*\)|use lib '$path2dir/lib\1|" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s|use lib.*|use lib '$path2dir/lib';|" $CURRENT/bin/nentferner.pl
$sedprog -i -e "s|use lib.*|use lib '$path2dir/lib';|g" $CURRENT/bin/translate.pl
$sedprog -i -e "s|use lib.*|use lib '$path2dir/lib';|g" $CURRENT/bin/oneSeq.pl

# paths to core_ortholog and blast_dir
echo "default path to blast_dir and core_orthologs: $path2dir"
$sedprog -i -e "s|\(my \$path = \).*|\1 '$path2dir';|g" $CURRENT/bin/hamstr.pl

###### CAN REMOVE THIS VAR $check_genewise in hamstr.pl ##########################
$sedprog -i -e 's/$check_genewise = [0,1];/$check_genewise = 1;/' $CURRENT/bin/hamstr.pl;
###############################################

### final check
echo "-------------------------------------"
echo "Final check..."
flag=0

echo "Perl modules"
for i in "${perlModules[@]}"; do
  msg=$((perl -e "use $i") 2>&1)
  if ! [[ -z ${msg} ]]; then
    echo "$i could not be installed"
    flag=1
  fi
done
echo "done!"

echo "Environment paths"
envPaths=(
  # "ONESEQDIR=$CURRENT"
  WISECONFIGDIR
)
for i in "${envPaths[@]}"; do
    if [ -z "$(grep $i ~/$bashFile)" ]; then
        echo "$i was not added into ~/$bashFile"
        flag=1
    fi
done
if [ "$fasta36" == "no" ]; then
    if [ -z "$(grep PATH=$CURRENT/bin/aligner/bin ~/$bashFile)" ]; then
        echo "$CURRENT/bin/aligner/bin was not added into ~/$bashFile"
        flag=1
    fi
fi
if [ -z "$(grep PATH=$CURRENT/bin:\$PATH ~/$bashFile)" ]; then
	echo "$CURRENT/bin was not added into ~/$bashFile"
fi

echo "done!"

if [ "$flag" == 1 ]; then
    echo "Some tools were not installed correctly or paths were not added into ~/$bashFile. Please run this setup again to try one more time!"
    exit
else
    echo "Generating symbolic link hamstr -> hamstr.pl"
    ln -s -f $CURRENT/bin/hamstr.pl $CURRENT/bin/hamstr
    echo "Sourcing bash profile file"
    source ~/$bashFile
    echo "All tests succeeded, HaMStR should be ready to run";
    $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/hamstr.pl;
    $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/oneSeq.pl;
    echo "Restart terminal and test your HaMStR with:"
    echo "perl bin/oneSeq.pl -seqFile=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global"
    echo "or"
    echo "perl bin/oneSeq.pl -h"
    echo "for more details."
fi
exit 1
