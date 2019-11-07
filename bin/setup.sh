#!/bin/bash

# check dependencies
echo "-------------------------------------"
echo "Checking dependencies..."

dependencies=(
  blastp
  genewise
  hmmsearch
  hmmbuild
  clustalw
  linsi
)

flag=0
for i in "${dependencies[@]}"; do
  msg="$(which $i)"
  if [ -z "$msg" ]; then
    echo "$i not found. Please install it to use HaMStR!"
    flag=1
  fi
done

if [ "$flag" == 1 ]; then exit 1; fi
echo "done!"

### prepare folder
echo "-------------------------------------"
echo "Preparing folders..."
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR/..

# create required folders
folders=(
  blast_dir
  core_orthologs
  genome_dir
  weight_dir
  taxonomy
  output
  tmp
)

for i in "${folders[@]}"; do
  if [ ! -d $i ]; then mkdir $i; fi
done

echo "done!"

### download data
echo "-------------------------------------"
echo "moving data into the right place"
echo "manipulate files in:"
CURRENT=$(pwd)
echo $CURRENT

if [[ $CURRENT == */HaMStR ]] || [[ $CURRENT == */hamstr ]]; then
    echo "Processing $CURRENT ..."
    echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar"
    wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/data_HaMStR.tar
    if [ ! -f $CURRENT/data_HaMStR.tar ]; then
        echo "File not found!"
    else
        CHECKSUM=$(cksum data_HaMStR.tar)
        echo "Checksum: $CHECKSUM"
        if [ "$CHECKSUM" == "557087663 1952579568 data_HaMStR.tar" ]; then
            echo "Extracting archive data_HaMStR.tar"
            tar xfv $CURRENT/data_HaMStR.tar
            echo "Archive data_HaMStR.tar extracted into $CURRENT"
            if [ ! -d $CURRENT/data_HaMStR ]; then
                echo "Directory $CURRENT/data_HaMStR not found!"
            else
                printf "\nMoving gene sets ...\n--------------------\n"
                #sleep 3
                rsync -rva data_HaMStR/genome_dir/* $CURRENT/genome_dir
                printf "\nMoving blast databases ...\n--------------------------\n"
                #sleep 3
                rsync -rva data_HaMStR/blast_dir/* $CURRENT/blast_dir
                printf "\nMoving annotations ...\n----------------------\n"
                #sleep 3
                rsync -rva data_HaMStR/weight_dir/* $CURRENT/weight_dir
                printf "\nMoving Taxonomy ...\n-------------------\n"
                #sleep 3
                rsync -rva data_HaMStR/taxonomy/* $CURRENT/taxonomy
                printf "\nMoving Pfam ...\n---------------\n"
                #sleep 3
                rsync -rva data_HaMStR/Pfam/* $CURRENT/bin/fas/Pfam
                printf "\nMoving SMART ...\n----------------\n"
                #sleep 3
                rsync -rva data_HaMStR/SMART/* $CURRENT/bin/fas/SMART
                printf "\nMoving CAST ...\n---------------\n"
                #sleep 3
                rsync -rva data_HaMStR/CAST/* $CURRENT/bin/fas/CAST
                printf "\nMoving COILS ...\n----------------\n"
                #sleep 3
                rsync -rva data_HaMStR/COILS2/* $CURRENT/bin/fas/COILS2
                printf "\nMoving SEG ...\n--------------\n"
                #sleep 3
                rsync -rva data_HaMStR/SEG/* $CURRENT/bin/fas/SEG
                printf "\nMoving SignalP ...\n------------------\n"
                #sleep 3
                rsync -rva data_HaMStR/SignalP/* $CURRENT/bin/fas/SignalP
                printf "\nMoving TMHMM ...\n----------------\n"
                #sleep 3
                rsync -rva data_HaMStR/TMHMM/* $CURRENT/bin/fas/TMHMM
                rsync -rva data_HaMStR/README* $CURRENT/
                printf "\nRemoving duplicated data. Please wait.\n------------------------------------\n"
                rm -rf $CURRENT/data_HaMStR
                printf "\nFinished. Data should be in place to run HaMStR.\n"
            fi
        else
            echo "Something went wrong with the download. Checksum does not match."
        fi
    fi

else
    echo "Please change into your HaMStR directory and run install_data.sh again."
    echo "Exiting."
    exit
fi
echo "done!"
