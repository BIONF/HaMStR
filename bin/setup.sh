#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

flag=0
### check grep, sed and wget availability
echo "-------------------------------------"
echo "Checking .bash_profile/.bashrc, grep, sed/gsed and wget availability..."
grepprog='grep'
sedprog='sed'
wgetprog='wget'
bashFile='.bashrc'
if [ "$sys" == "Darwin" ]; then
    sedprog='gsed'
    bashFile='.bash_profile'
fi

# NOTE: install only available for Linux!
if [ -z "$(which $sedprog)" ]; then
    if [ "$sys" == "Darwin" ]; then
        echo "$sedprog not found. Please install it first (e.g. using brew)!"
        flag=1
    fi
fi

if [ -z "$(which $grepprog)" ]; then
    if [ "$sys" == "Darwin" ]; then
        echo "$grepprog not found. Please install it first (e.g. using brew)!"
        flag=1
    fi
fi

if [ -z "$(which $wgetprog)" ]; then
    if [ "$sys" == "Darwin" ]; then
        echo "$wgetprog not found. Please install it first (e.g. using brew)!"
        flag=1
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

if [ -z "$(which R)" ]; then
    echo "R"
    conda install -y r
fi

# if [[ -z $(conda list | grep "pkg-config ") ]]; then
#     echo "pkg-config"
#     conda install -y pkg-config
# fi
#
# if [[ -z $(conda list | grep "perl-bioperl ") ]]; then
#     echo "perl-bioperl"
#     conda install -y -c bioconda perl-bioperl
#     conda install -y -c bioconda perl-bioperl-core
#     conda install -y -c bioconda perl-bioperl-run
# fi

dependencies=(
  blastp # blast
  blastall # blast-legacy
  genewise # wise2
  hmmsearch # hmmer (for both hmmsearch and hmmbuild)
  clustalw
  mafft # for linsi
  muscle
)

# for i in "${dependencies[@]}"; do
#   if [ -z "$(which $i)" ]; then
#     echo $i
#     tool=$i
#     if [ "$tool" = "blastp" ]; then
#       conda install -y -c bioconda blast
#     elif [ "$tool" = "blastall" ]; then
#       conda install -y -c bioconda blast-legacy
#     elif [ "$tool" = "hmmsearch" ]; then
#       conda install -y -c bioconda hmmer
#     elif [ "$tool" = "genewise" ]; then
#       conda install -y -c bioconda wise2
#       wisePath=$(which "genewise")
#       if [ -z "$(grep WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
#           echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
#       fi
#     else
#       conda install -y -c bioconda $i
#     fi
#   fi
# done

for i in "${dependencies[@]}"; do
  if [ -z "$(which $i)" ]; then
    echo "$i not found. Please install it to use HaMStR!"
    flag=1
  fi
done
if [ "$flag" == 1 ]; then exit 1; fi

perlModules=(
  DBI
  DB_File
  File::Copy
  File::Path
  File::Basename
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

for i in "${perlModules[@]}"; do
  msg=$((perldoc -l $i) 2>&1)
  if [[ "$(echo $msg)" == *"No documentation"* ]]; then
    cpanm ${i} --force
  fi
done

echo "done!"

# ### prepare folder
# echo "-------------------------------------"
# echo "Preparing folders..."
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# cd $DIR/..
# CURRENT=$(pwd)
#
# # create required folders
# folders=(
#   blast_dir
#   core_orthologs
#   genome_dir
#   weight_dir
#   taxonomy
#   output
#   "bin/fas/CAST"
#   "bin/fas/COILS2"
#   "bin/fas/Pfam"
#   "bin/fas/Pfam/Pfam-hmms"
#   "bin/fas/Pfam/output_files"
#   "bin/fas/SEG"
#   "bin/fas/SignalP"
#   "bin/fas/SMARublastT"
#   "bin/fas/TMHMM"
#   "bin/aligner"
# )
#
# for i in "${folders[@]}"; do
#   echo "$i"
#   if [ ! -d $i ]; then mkdir $i; fi
# done
# echo "done!"
#
# ### download tools
# echo "-------------------------------------"
# echo "Downloading and installing annotation tools/databases:"
#
# fasta36="yes"
# if [ -z "$(which fasta36)" ]; then
#   fasta36="no"
#   echo "fasta-36"
#   fasta36v="fasta-36.3.8h"
#   wget "http://faculty.virginia.edu/wrpearson/fasta/fasta36/${fasta36v}.tar.gz"
#   tar xfv $fasta36v.tar.gz
#   rm "${fasta36v}.tar.gz"
#   mv $fasta36v bin/aligner
#   cd "bin/aligner/$fasta36v/src"
#   if [ $sys=="Linux" ]; then
#     make -f ../make/Makefile.linux64_sse2 all
#   elif [ $sys=="Darwin" ]; then
#     make -f ../make/Makefile.os_x86_64 all
#   fi
#   fastaPath=$(cd -- ../bin && pwd)
#   if [ -z "$(grep PATH=${fastaPath} ~/$bashFile)" ]; then
#       echo "export PATH=${fastaPath}:\$PATH" >> ~/$bashFile
#   fi
# fi
# cd $CURRENT
#
# seg="yes"
# if [ -z "$(which seg)" ]; then
#   seg="no"
#   echo "SEG"
#   cd "bin/fas/SEG"
#   wget -r -l 2 -np ftp://ftp.ncbi.nih.gov/pub/seg/seg
#   mv ftp.ncbi.nih.gov/pub/seg/seg/* $(pwd)
#   rm -rf ftp.ncbi.nih.gov
#   rm -rf archive
#   make
#   segPath=$(pwd)
#   if [ -z "$(grep PATH=${segPath} ~/$bashFile)" ]; then
#       echo "export PATH=${segPath}:\$PATH" >> ~/$bashFile
#   fi
# fi
# cd $CURRENT
#
# cd "bin/fas/Pfam/Pfam-hmms"
# if ! [ -f Pfam-A.hmm ]; then
#   echo "pfam-A.hmm"
#   wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release//Pfam-A.hmm.gz
#   wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
#   gunzip Pfam-A.hmm.gz
#   gunzip Pfam-A.hmm.dat.gz
#   hmmpress Pfam-A.hmm
#   mv $CURRENT/bin/pfam_scan.pl $CURRENT/bin/fas/Pfam/
# fi
# cd $CURRENT
#
# ### download data
# echo "-------------------------------------"
# echo "Getting pre-calculated data"
#
# if ! [ "$(ls -A $CURRENT/taxonomy)" ]; then
#     # if [[ $CURRENT == */HaMStR ]] || [[ $CURRENT == */hamstr ]]; then
#       echo "Processing $CURRENT ..."
#       echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
#       wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar
#       if [ ! -f $CURRENT/data_HaMStR.tar ]; then
#         echo "File data_HaMStR.tar not found! Please try again!"
#       else
#         CHECKSUM=$(cksum data_HaMStR.tar)
#         echo "Checksum: $CHECKSUM"
#         if [ "$CHECKSUM" == "557087663 1952579568 data_HaMStR.tar" ]; then
#           echo "Extracting archive data_HaMStR.tar"
#           tar xfv $CURRENT/data_HaMStR.tar
#           rm $CURRENT/data_HaMStR.tar
#           echo "Archive data_HaMStR.tar extracted into $CURRENT"
#           if [ ! -d $CURRENT/data_HaMStR ]; then
#             echo "Directory $CURRENT/data_HaMStR not found!"
#           else
#             printf "\nMoving gene sets ...\n--------------------\n"
#             rsync -rva data_HaMStR/genome_dir/* $CURRENT/genome_dir
#             printf "\nMoving blast databases ...\n--------------------------\n"
#             rsync -rva data_HaMStR/blast_dir/* $CURRENT/blast_dir
#             printf "\nMoving annotations ...\n----------------------\n"
#             rsync -rva data_HaMStR/weight_dir/* $CURRENT/weight_dir
#             printf "\nMoving Taxonomy ...\n-------------------\n"
#             rsync -rva data_HaMStR/taxonomy/* $CURRENT/taxonomy
#             # printf "\nMoving Pfam ...\n---------------\n"
#             # rsync -rva data_HaMStR/Pfam/* $CURRENT/bin/fas/Pfam
#             printf "\nMoving SMART ...\n----------------\n"
#             rsync -rva data_HaMStR/SMART/* $CURRENT/bin/fas/SMART
#             printf "\nMoving CAST ...\n---------------\n"
#             rsync -rva data_HaMStR/CAST/* $CURRENT/bin/fas/CAST
#             printf "\nMoving COILS ...\n----------------\n"
#             rsync -rva data_HaMStR/COILS2/* $CURRENT/bin/fas/COILS2
#             # printf "\nMoving SEG ...\n--------------\n"echo "export ONESEQDIR=${CURRENT}" >> ~/$bashFile
#             # rsync -rva data_HaMStR/SEG/* $CURRENT/bin/fas/SEG
#             printf "\nMoving SignalP ...\n------------------\n"
#             rsync -rva data_HaMStR/SignalP/* $CURRENT/bin/fas/SignalP
#             printf "\nMoving TMHMM ...\n----------------\n"
#             rsync -rva data_HaMStR/TMHMM/* $CURRENT/bin/fas/TMHMM
#             rsync -rva data_HaMStR/README* $CURRENT/
#             printf "\nRemoving duplicated data. Please wait.\n------------------------------------\n"
#             rm -rf $CURRENT/data_HaMStR
#             printf "\nDone! Data should be in place to run HaMStR.\n"
#           fi
#         else
#           echo "Something went wrong with the download. Checksum does not match."
#         fi
#       fi
#     # else
#     #   echo "Please change into your HaMStR directory and run install_data.sh again."
#     #   echo "Exiting."
#     #   exit
#     # fi
# fi
#
# ### add paths to bash profile file
# echo "-------------------------------------"
# echo "Adding paths to ~/$bashFile"
#
# if [ -z "$(grep ONESEQDIR=$CURRENT ~/$bashFile)" ]; then
#   echo "export ONESEQDIR=${CURRENT}" >> ~/$bashFile
# fi
#
# wisePath=$(which "genewise")
# if [ -z "$(grep WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
#     echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
# fi
# echo "done!"
#
# ### adapt paths in hamstr scripts
# echo "-------------------------------------"
# echo "Adapting paths in hamstr scripts"
# # update the sed and grep commands
# $sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/hamstr.pl;
# $sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/hamstr.pl;
# $sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/oneSeq.pl;
# $sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/oneSeq.pl;
#
# # localize the perl installation
# path2perl=`which perl`
# echo "path to perl: $path2perl"
# $sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/hamstr.pl;
# $sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/nentferner.pl;
# $sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/translate.pl;
# $sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/oneSeq.pl;
#
# # get lib path
# path2dir=$CURRENT
# echo "path to lib: $path2dir/lib"
# $sedprog -i -e "s|use lib.*lib\(.*\)|use lib '$path2dir/lib\1|" $CURRENT/bin/hamstr.pl
# $sedprog -i -e "s|use lib.*|use lib '$path2dir/lib';|" $CURRENT/bin/nentferner.pl
# $sedprog -i -e "s|use lib.*|use lib '$path2dir/lib';|g" $CURRENT/bin/translate.pl
# $sedprog -i -e "s|use lib.*|use lib '$path2dir/lib';|g" $CURRENT/bin/oneSeq.pl
#
# # paths to core_ortholog and blast_dir
# echo "default path to blast_dir and core_orthologs: $path2dir"
# $sedprog -i -e "s|\(my \$path = \).*|\1 '$path2dir';|g" $CURRENT/bin/hamstr.pl
#
# ###### CAN REMOVE THIS VAR $check_genewise in hamstr.pl ##########################
# $sedprog -i -e 's/$check_genewise = [0,1];/$check_genewise = 1;/' $CURRENT/bin/hamstr.pl;
# ###############################################
#
# ### final check
# echo "-------------------------------------"
# echo "Final check..."
# flag=0
#
# echo "Conda packages"
# condaPkgs=(
#   perl-bioperl
#   perl-bioperl-core
#   blast
#   blast-legacy
#   hmmer
#   wise2
#   clustalw
#   mafft
#   muscle
# )
# for i in "${condaPkgs[@]}"; do
#     if [[ -z $(conda list | grep "$i ") ]]; then
#         echo "$i could not be installed"
#         flag=1
#     fi
# done
# echo "done!"
#
# echo "Perl modules"
# for i in "${perlModules[@]}"; do
#   msg=$((perl -e "use $i") 2>&1)
#   if ! [[ -z ${msg} ]]; then
#     echo "$i could not be installed"
#     flag=1
#   fi
# done
# echo "done!"
#
# echo "Environment paths"
# envPaths=(
#   "ONESEQDIR=$CURRENT"
#   WISECONFIGDIR
# )
# for i in "${envPaths[@]}"; do
#     if [ -z "$(grep $i ~/$bashFile)" ]; then
#         echo "$i was not added into ~/$bashFile"
#         flag=1
#     fi
# done
# if [ "$fasta36" == "no" ]; then
#     if [ -z "$(grep PATH=$CURRENT/bin/aligner/fasta-36 ~/$bashFile)" ]; then
#         echo "$CURRENT/bin/aligner/fasta-36 was not added into ~/$bashFile"
#         flag=1
#     fi
# fi
# if [ "$seg" == "no" ]; then
#     if [ -z "$(grep PATH=$CURRENT/bin/fas/SEG ~/$bashFile)" ]; then
#         echo "$CURRENT/bin/fas/SEG was not added into ~/$bashFile"
#         flag=1
#     fi
# fi
# echo "done!"
#
# if [ "$flag" == 1 ]; then
#     echo "Some tools were not installed correctly or paths were not added into ~/$bashFile. Please check again!"
#     exit
# else
#     echo "Generating symbolic link hamstr -> hamstr.pl"
#     ln -s -f $CURRENT/bin/hamstr.pl $CURRENT/bin/hamstr
#     echo "Sourcing bash profile file"
#     source ~/$bashFile
#     echo "All tests succeeded, HaMStR should be ready to run";
#     $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/hamstr.pl;
#     $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/oneSeq.pl;
#     echo "Restart terminal and test your HaMStR with:"
#     echo "perl bin/oneSeq.pl -sequence_file=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global"
#     echo "or"
#     echo "perl bin/oneSeq.pl -h"
#     echo "for more details."
# fi
# exit 1
