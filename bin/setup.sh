#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

flag=0
root=1
fas=1

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

if [ "$EUID" -ne 0 ]; then
    echo "You are not running this setup as root."
    read -p "Press enter to continue, but some missing tools/libraries will not be installed!"
    root=0
fi

### check grep, sed and wget availability
echo "-------------------------------------"
echo "Checking .bash_profile/.bashrc, grep, sed/gsed and wget availability..."
grepprog='grep'
sedprog='sed'
wgetprog='wget'
bashFile='.bashrc'
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
    exit
fi

if [ -z "$(which $grepprog)" ]; then
    echo -e "\033[31m$grepprog not found!\033[0m"
    exit
fi

if [ -z "$(which $wgetprog)" ]; then
    echo -e "\033[31m$wgetprog not found!\033[0m"
    exit
fi

if ! [ -f ~/$bashFile ]; then
    touch ~/$bashFile
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

### install dependencies
if [ $root == 1 ]; then
    sudo bash $CURRENT/install_lib.sh
fi

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
    rm *.dmp
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
                annoFAS --fasta test.fa --path $CURRENT --name q --prepare 1 --annoPath $CURRENT/bin/fas
            fi
        else
            pip install $CURRENT/bin/fas --user
            if [ -z "$(grep \$HOME/.local/bin:\$PATH ~/$bashFile)" ]; then
                echo "export PATH=\$HOME/.local/bin:\$PATH" >> ~/$bashFile
            fi
            # change path to annoFAS.py and greeyFAS.py in oneSeq.pl (to not require for restarting the terminal)
            annoprog="python \$path\/bin\/fas\/greedyFAS\/annoFAS.py"
            $sedprog -i -e "s/\(my \$annotation_prog = \).*/\1\"$annoprog\";/" $CURRENT/bin/oneSeq.pl
            fasprog="python \$path\/bin\/fas\/greedyFAS\/greedyFAS.py"
            $sedprog -i -e "s/\(my \$fas_prog = \).*/\1\"$fasprog\";/" $CURRENT/bin/oneSeq.pl
            # get FAS annotation tools and pre-calculated data
            python $CURRENT/bin/fas/greedyFAS/annoFAS.py --fasta $CURRENT/data/infile.fa --path $CURRENT --name q --prepare 1 --annoPath $CURRENT/bin/fas
        fi
    else
        fasPath="$(pip show greedyFAS | grep Location | sed 's/Location: //')"
        annoFile="$fasPath/greedyFAS/annoFAS.pl"
        tmp="$(grep "my \$config" $annoFile | sed 's/my \$config = //' | sed 's/;//')"
        if [ $tmp == "1" ]; then
            annoPath="$(grep "my \$annotationPath" $annoFile | sed 's/my \$annotationPath = "//' | sed 's/";//')"
            echo "$annoPath"
            if ! [ -f "$annoPath/Pfam/Pfam-hmms/Pfam-A.hmm" ]; then
                annoFAS --fasta $CURRENT/data/infile.fa --path $CURRENT --name q --prepare 1 --annoPath $annoPath
            fi
        else
            annoFAS --fasta $CURRENT/data/infile.fa --path $CURRENT --name q --prepare 1 --annoPath $CURRENT/bin/fas
        fi
    fi

    cd $CURRENT
    echo "done!"
fi

### download data
echo "-------------------------------------"
echo "Getting pre-calculated data"

data_HaMStR_file="data_HaMStR2018b.tar.gz"
checkSumData="979235026 675298057 $data_HaMStR_file"

if ! [ "$(ls -A $CURRENT/genome_dir)" ]; then
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
  blastn
)

for i in "${dependencies[@]}"; do
  if [ -z "$(which $i)" ]; then
    echo -e "\t\033[31mWARNING $i not found!\033[0m"
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
echo "done!"

echo "Perl modules"
for i in "${perlModules[@]}"; do
  msg=$((perl -e "use $i") 2>&1)
  if ! [[ -z ${msg} ]]; then
    echo -e "\t\033[31mWARNING $i could not be installed\033[0m"
    flag=1
  fi
done
echo "done!"

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
echo "done!"

if [ "$flag" == 1 ]; then
    echo "Some tools/libraries could not be found or paths were not added into ~/$bashFile."
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
    echo "All tests succeeded, HaMStR should be ready to run"
    $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/hamstr.pl
    $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/oneSeq.pl
    echo "Test your HaMStR with:"
    echo -e "\033[1mperl bin/oneSeq.pl -seqFile=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global\033[0m"
    echo "or"
    echo -e "\033[1mperl bin/oneSeq.pl -h\033[0m"
    echo "for more details."
fi
exit 1
