#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

CURRENT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

flag=0
root=0
fas=1
homedir="$(echo $HOME)"

while getopts ":f" opt; do
    case ${opt} in
        f )
          echo "NO FAS!"
          fas=0
          ;;
        \? )
          echo "wrong option"
          exit 1
          ;;
    esac
done

if [ "$EUID" -eq 0 ]; then
    echo "Please DO NOT run this script as root!"
    # read -p "Press enter to continue, but some missing tools/libraries will not be installed!"
    exit
    # root=0
else
    read -p "Do you have sudo password? [y/n]" -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        root=1
    fi
fi

### install dependencies
if [ "$sys" == "Darwin" ]; then
    bash $CURRENT/install_lib.sh
else
    if [ $root == 1 ]; then
        echo "Enter sudo password to install required libraries..."
        sudo bash $CURRENT/install_lib.sh
    fi
fi

### check grep, sed and wget availability
echo "-------------------------------------"
echo "Checking .bash_profile/.bashrc, grep, sed/gsed and wget availability..."
grepprog='grep'
sedprog='sed'
wgetprog='wget'
bashFile='.bashrc'
rprofile='.Rprofile'

if [ "$sys" == "Darwin" ]; then
    sedprog='gsed'
	grepprog='ggrep'
	shell=$(echo $SHELL)
	if [ $shell == "/bin/zsh" ]; then
    	bashFile='.zshrc'
	else
		bashFile='.bash_profile'
	fi
fi

if [ -z "$(which $sedprog)" ]; then
    echo -e "\033[31m$sedprog not found!\033[0m"
	if [ $root == 0 ]; then
		echo "Please run $CURRENT/install_lib.sh first!"
    	exit
	fi
fi

if [ -z "$(which $grepprog)" ]; then
    echo -e "\033[31m$grepprog not found!\033[0m"
	if [ $root == 0 ]; then
		echo "Please run $CURRENT/install_lib.sh first!"
    	exit
	fi
fi

if [ -z "$(which $wgetprog)" ]; then
    echo -e "\033[31m$wgetprog not found!\033[0m"
	if [ $root == 0 ]; then
		echo "Please run $CURRENT/install_lib.sh first!"
    	exit
	fi
fi

if ! [ -f ~/$bashFile ]; then
    touch ~/$bashFile
fi
if ! [ -f ~/$rprofile ]; then
    touch ~/$rprofile
fi
echo "done!"

### prepare folders
echo "-------------------------------------"
echo "Preparing folders..."

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
    if [ ! -d "$CURRENT/$i" ]; then mkdir "$CURRENT/$i"; fi
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
    	  tar xf $fasta36v.tar.gz
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
		echo -e "\033[31mfasta36 tool could not be found in $CURRENT/bin/aligner/. Please check again!\033[0m"
		exit
	fi
fi

cd "taxonomy"
if ! [ -f "nodes" ]; then
    wget "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tar xf taxdump.tar.gz
    rm taxdump.tar.gz
    echo "Taxonomy database indexing. It can take a while, please wait..."
    perl $CURRENT/bin/indexTaxonomy.pl $CURRENT/taxonomy
    rm citations.dmp
    rm delnodes.dmp
    rm division.dmp
    rm gencode.dmp
    rm merged.dmp
    rm gc.prt
    rm readme.txt
fi
cd $CURRENT
if ! [ -f "$CURRENT/taxonomy/nodes" ]; then
	echo -e "\033[31mError while indexing NCBI taxonomy database! Please check $CURRENT/taxonomy/ folder and run this setup again!\033[0m"
	exit
fi

if [ $fas == 1 ]; then
    cd "bin"
    if [ -z "$(which greedyFAS)" ]; then
        echo "FAS"
        if ! [ -f "fas/setup.py" ]; then
            wget https://github.com/BIONF/FAS/archive/master.tar.gz
            tar xf master.tar.gz
            mv FAS-master fas
            rm master.tar.gz
        fi
        if [ $root == 1 ]; then
            pip install $CURRENT/bin/fas
            if [ -z "$(which annoFAS)" ]; then
                echo "Installation of FAS failed! Please try again!"
                exit
            else
                annoFAS --fasta test.fa --path $CURRENT --name q --prepare --annoPath $CURRENT/bin/fas
            fi
        else
            pip install $CURRENT/bin/fas --user
            if [ -z "$(grep \$HOME/.local/bin:\$PATH ~/$bashFile)" ]; then
                echo "export PATH=\$HOME/.local/bin:\$PATH" >> ~/$bashFile
            fi
            if [ -z "$(grep $homedir/.local/bin ~/$rprofile)" ]; then
                echo "Sys.setenv(PATH = paste(\"$homedir/.local/bin\", Sys.getenv(\"PATH\"), sep=\":\"))" >> ~/$rprofile
            fi
            # change path to annoFAS.py and greeyFAS.py in oneSeq.pl (to not require for restarting the terminal)
            annoprog="python \$path\/bin\/fas\/greedyFAS\/annoFAS.py"
            $sedprog -i -e "s/\(my \$annotation_prog = \).*/\1\"$annoprog\";/" $CURRENT/bin/oneSeq.pl
            fasprog="python \$path\/bin\/fas\/greedyFAS\/greedyFAS.py"
            $sedprog -i -e "s/\(my \$fas_prog = \).*/\1\"$fasprog\";/" $CURRENT/bin/oneSeq.pl
            # get FAS annotation tools and pre-calculated data
            python $CURRENT/bin/fas/greedyFAS/annoFAS.py --fasta $CURRENT/data/infile.fa --path $CURRENT --name q --prepare --annoPath $CURRENT/bin/fas
        fi
    else
        fasPath="$(pip show greedyFAS | grep Location | sed 's/Location: //')"
        annoFile="$fasPath/greedyFAS/annoFAS.pl"
        tmp="$(grep "my \$config" $annoFile | sed 's/my \$config = //' | sed 's/;//')"
        if [ $tmp == "1" ]; then
            annoPath="$(grep "my \$annotationPath" $annoFile | sed 's/my \$annotationPath = "//' | sed 's/";//')"
            if ! [ -f "$annoPath/Pfam/Pfam-hmms/Pfam-A.hmm" ]; then
                annoFAS --fasta $CURRENT/data/infile.fa --path $CURRENT --name q --prepare --annoPath $annoPath
            fi
        else
            annoFAS --fasta $CURRENT/data/infile.fa --path $CURRENT --name q --prepare --annoPath $CURRENT/bin/fas
        fi
    fi

    cd $CURRENT
    echo "done!"
fi

### download data
data_HaMStR_file="data_HaMStR-2019.tar.gz"
checkSumData="1303809705 685885017 $data_HaMStR_file"

if ! [ "$(ls -A $CURRENT/genome_dir)" ]; then
    echo "-------------------------------------"
    echo "Getting pre-calculated data"

	echo "Processing $CURRENT ..."
	if [ ! -f $CURRENT/$data_HaMStR_file ]; then
		echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_HaMStR_file"
		wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_HaMStR_file
	else
		CHECKSUM=$(cksum $data_HaMStR_file)
		echo "Checksum: $CHECKSUM"
		if ! [ "$CHECKSUM" == "$checkSumData" ]; then
    		  rm $CURRENT/$data_HaMStR_file
    		  echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_HaMStR_file"
      		  wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_HaMStR_file
    	fi
    fi

	if [ ! -f $CURRENT/$data_HaMStR_file ]; then
        echo "File $data_HaMStR_file not found! Please try to download again from"
        echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
        exit
	fi

	CHECKSUM=$(cksum $data_HaMStR_file)
	if [ "$CHECKSUM" == "$checkSumData" ]; then
        echo "Extracting archive $data_HaMStR_file..."
        tar xf $CURRENT/$data_HaMStR_file
        rm $CURRENT/$data_HaMStR_file
        for i in $(ls "$CURRENT/genome_dir"); do rm "$CURRENT/genome_dir/$i/$i.fa.mod"; done

        if [ "$(ls -A $CURRENT/blast_dir)" ]; then
            echo "Data should be in place to run HaMStR."
        else
            echo -e "\033[31mSomething went wrong with the download. Data folders are empty.\033[0m"
            echo "Please try to download again from"
            echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_HaMStR_file"
            echo "Or contact us if you think this is our issue!"
            exit
        fi
	else
        echo -e "\033[31mSomething went wrong with the download. Checksum does not match.\033[0m"
        echo "Please try to download again from"
        echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_HaMStR_file"
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

echo "Adding paths to ~/$rprofile"
if [ -z "$(grep $CURRENT/bin ~/$rprofile)" ]; then
    echo "Sys.setenv(PATH = paste(\"$CURRENT/bin\", Sys.getenv(\"PATH\"), sep=\":\"))" >> ~/$rprofile
fi
echo "done!"

### adapt paths in hamstr scripts
echo "-------------------------------------"
echo "Adapting paths in hamstr scripts"
# update the sed and grep commands
$sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/oneSeq.pl
$sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/oneSeq.pl

# localize the perl installation
path2perl=`which perl`
echo "path to perl: $path2perl"
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/nentferner.pl
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/translate.pl
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/oneSeq.pl

echo "done!"

### final check
echo "-------------------------------------"
echo "Final check..."
flag=0

echo "Tools"
dependencies=(
  genewise
  hmmsearch
  hmmscan
  hmmbuild
  mafft
  muscle
  clustalw
  blastp
)

for i in "${dependencies[@]}"; do
    tool=$i
    if [ $tool == "clustalw" ]; then
        if [ "$sys" == "Darwin" ]; then
            tool="clustalw2"
        fi
    fi
    if [ -z "$(which $tool)" ]; then
        echo -e "\t\033[31mWARNING $tool not found!\033[0m"
        flag=1
    fi
done

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

echo "Perl modules"
for i in "${perlModules[@]}"; do
  msg=$((perl -e "use $i") 2>&1)
  if ! [[ -z ${msg} ]]; then
    echo -e "\t\033[31mWARNING $i could not be installed\033[0m"
    flag=1
  fi
done

echo "Environment paths"
envPaths=(
  WISECONFIGDIR
)
for i in "${envPaths[@]}"; do
    if [ -z "$(grep $i ~/$bashFile)" ]; then
        echo -e "\t\033[31mWARNING $i was not added into ~/$bashFile\033[0m"
        flag=1
    fi
done
if [ "$fasta36" == "no" ]; then
    if [ -z "$(grep PATH=$CURRENT/bin/aligner/bin ~/$bashFile)" ]; then
        echo -e "\t\033[31mWARNING $CURRENT/bin/aligner/bin was not added into ~/$bashFile\033[0m"
        flag=1
    fi
fi
if [ -z "$(grep PATH=$CURRENT/bin:\$PATH ~/$bashFile)" ]; then
	echo -e "\t\033[31mWARNING $CURRENT/bin was not added into ~/$bashFile\033[0m"
    flag=1
fi
if [ -z "$(grep $CURRENT/bin ~/$rprofile)" ]; then
	echo -e "\t\033[31mWARNING $CURRENT/bin was not added into ~/$rprofile\033[0m"
    flag=1
fi

if [ "$flag" == 1 ]; then
    echo "Some tools/libraries could not be found or paths were not added into ~/$bashFile or ~/$rprofile."
    echo "Please install the missing dependencies using $CURRENT/install_lib.sh script (or ask your admin if you don't have root privileges)."
    echo "Then run this setup again to try one more time!"
    exit
else
    echo "Generating symbolic links"
    ln -s -f $CURRENT/bin/hamstr.pl $CURRENT/bin/hamstr
    ln -s -f $CURRENT/bin/oneSeq.pl $CURRENT/bin/oneSeq
    echo "Sourcing bash profile file"
    source ~/$bashFile
    echo "-------------------------------------"
    $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/hamstr.pl
    $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/oneSeq.pl
    echo "All tests succeeded, HaMStR should be ready to run. You can test it with:"
    echo -e "\033[1moneSeq -seqFile=infile.fa -seqName=test -refspec=HUMAN@9606@3 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global -cpu=4\033[0m"
    echo "Output files with prefix \"test\" will be found at your current working directory!"
    echo "For more details, use"
    echo -e "\033[1moneSeq -h\033[0m"
fi
exit 1
