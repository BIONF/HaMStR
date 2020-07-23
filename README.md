# HaMStR-OneSeq
[![conda-install](https://anaconda.org/bionf/hamstr/badges/installer/conda.svg)](https://anaconda.org/bionf/hamstr)
[![conda-version](https://anaconda.org/bionf/hamstr/badges/version.svg)](https://anaconda.org/bionf/hamstr)
[![GPLv3-license](https://anaconda.org/bionf/hamstr/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.de.html)

# Table of Contents
* [How to install](#how-to-install)
     * [0. Basic system tools requirement](#0-basic-system-tools-requirement)
     * [1. Dependencies](#1-dependencies)
     * [2a. Install using Anaconda](#1a-install-using-anaconda)
     * [2b. Install in Ubuntu/MacOS](#1b-install-in-ubuntumacos)
* [Usage](#usage)
* [HaMStR data set](#hamstr-data-set)
* [Bugs](#bugs)
* [How to cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

## 0. Basic system tools requirement
You need to have `wget`, `grep` and `sed` (or `gsed` for **MacOS**) to install HaMStR. So please install them if they are missing. For MacOS users, we recommend using [Homebrew](https://brew.sh) to install those command line tools.
To use [FAS tool](https://github.com/BIONF/FAS) (a dependency of HaMStR), you also need [Python 3](https://www.python.org/downloads/).

## 1. Dependencies

*HaMStR-oneSeq* has some dependencies, that either will be automatically installed via the setup script, or must be installed by your system admin if you don't have the root privileges. In [our wiki](https://github.com/BIONF/HaMStR/wiki/Dependencies) you will find the full list of HaMStR-oneSeq's dependencies for Ubuntu system as well as the alternatives for MacOS. 

In Ubuntu, you can install those system and bioinformatics tools/libraries using `apt-get` tool
```
sudo apt-get update -y
sudo apt-get install tool_name -y
```
In MacOS, we suggest using [Homebrew](https://brew.sh) as a replacement for `apt-get`. After having Homebrew, you can install tools/libraries by using the command
```
brew install tool_name
```
In both operation systems, you can install Perl modules using `cpanm`.
```
# first, install cpanm
curl -L http://cpanmin.us | perl - --sudo App::cpanminus
# then, install perl module using cpanm
sudo cpanm perl_module_name
```

If you do not have root privileges, ask your admin to install these dependencies using the `install_lib.sh` script.

```
cd HaMStR
sudo ./install_lib.sh
```

## 2a. Install using Anaconda

Follow [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install conda (anaconda or miniconda) to your system.

Add additional channels [bioconda](https://bioconda.github.io/) and [conda-forge](https://conda-forge.org/):
```
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create and activate a conda environment for HaMStR
```
conda create --name hamstr -y
conda activate hamstr
```

Install HaMStR
```
conda install -c BIONF hamstr
setup_hamstr
```
HaMStR will be installed under the subfolder **HaMStR** in side your current working directory.
After the setup run successfully, you can start using HaMStR (in some cases you should restart the terminal).

## 2b. Install in Ubuntu/MacOS

Get HaMStR source code from GitHub
```
git clone --depth=1 https://github.com/BIONF/HaMStR
```

Run `setup.sh` script in the HaMStR folder to install HaMStR and its dependencies
```
cd HaMStR
./setup.sh
```
*You should have the sudo password ready, otherwise some missing dependencies cannot be installed. See [dependency list](#dependencies) for more info. If you do not have root privileges, ask your admin to install those dependencies using `install_lib.sh` script.*

After the setup run successfully, you can start using HaMStR (in some cases you should restart the terminal).

*For debugging the installation, please create a log file by running the setup as e.g. `bin/setup.sh | tee log.txt` for Linux/MacOS or `setup_hamstr | tee log.txt` for Anaconda and send us that log file, so that we can trouble shoot the issues. Most of the problems can be solved by just re-running the setup.*

# Usage
HaMStR will run smoothly with the provided sample input file in 'HaMStR/data/infile.fa' if everything is set correctly.

```
oneSeq -seqFile=infile.fa -seqName=test -refspec=HUMAN@9606@3 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -cpu=8
```
The output files with the prefix `test` will be saved at your current working directory.
You can have an overview about the available options with the command
```
oneSeq -h
```

*If you get the error message that `oneSeq command not found`, you should restart the terminal, or replace `oneSeq` by `perl bin/oneSeq`*

Please find more information in [our wiki](https://github.com/BIONF/HaMStR/wiki) to learn about the [input and outputs files](https://github.com/BIONF/HaMStR/wiki/Input-and-Output-Files) of *HaMStR-oneSeq*. 

# HaMStR data set

Within the data package (https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) we provide a set of 78 reference taxa. They can be automatically downloaded during the setup. This data comes "ready to use" with the HaMStR-OneSeq framework. Species data must be present in the three directories listed below:

* genome_dir (Contains sub-directories for proteome fasta files for each species)
* blast_dir (Contains sub-directories for BLAST databases made with `makeblastdb` out of your proteomes)
* weight_dir (Contains feature annotation files for each proteome)

For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version])

HaMStR-oneSeq is not limited to those 78 taxa. If needed the user can manually add further gene sets (multifasta format) using provided python scripts.

## Adding a new gene set into HaMStR
For adding **one gene set**, please use the `bin/addTaxonHamstr.py` script:
```
python3 bin/addTaxonHamstr.py -f newTaxon.fa -i tax_id -o /path/to/HaMStR [-n abbr_tax_name] [-c] [-v protein_version] [-a]
```

in which, the first 3 arguments are required including `newTaxon.fa` is the gene set that need to be added, `tax_id` is its NCBI taxonomy ID, `/path/to/HaMStR` is where the sub-directories can be found (*genome_dir*, *blast_dir* and *weight_dir*). Other arguments are optional, which are `-n` for specify your own taxon name (if not given, an abbriviate name will be suggested based on the NCBI taxon name of the input `tax_id`), `-c` for calculating the BLAST DB (only needed if you need to include your new taxon into the list of taxa for compilating the core set), `-v` for identifying the genome/proteome version (default will be 1), and `-a` for turning off the annotation step (*not recommended*).

## Adding a list of gene sets into HaMStR
For adding **more than one gene set**, please use the `bin/addTaxaHamstr.py` script:
```
python3 bin/addTaxaHamstr.py -i /path/to/newtaxa/fasta -m mapping_file -o /path/to/HaMStR [-c]
```
in which, `/path/to/taxa/fasta` is a folder where the FASTA files of all new taxa can be found. `mapping_file` is a tab-delimited text file, where you provide the taxonomy IDs that stick with the FASTA files:

```
#filename	tax_id	abbr_tax_name	version
filename1	12345678
filename2	9606
filename3	4932	my_fungi
...
```

The header line (started with #) is a Must. The values of the last 2 columns (abbr. taxon name and genome version) are, however, optional. If you want to specify a new version for a genome, you need to define also the abbr. taxon name, so that the genome version is always at the 4th column in the mapping file.

_**NOTE:** After adding new taxa into *HaMStR-oneSeq*, you should [check for the validity of the new data](https://github.com/BIONF/HaMStR/wiki/Check-data-validity) before running HaMStR._

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/HaMStR/issues/new) or be in touch via email.

# How to cite
Ebersberger, I., Strauss, S. & von Haeseler, A. HaMStR: Profile hidden markov model based search for orthologs in ESTs. BMC Evol Biol 9, 157 (2009), [doi:10.1186/1471-2148-9-157](https://doi.org/10.1186/1471-2148-9-157)

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Vinh Tran](https://github.com/trvinh)
- [Holger Bergmann](https://github.com/holgerbgm)

# Contact
For further support or bug reports please contact: ebersberger@bio.uni-frankfurt.de
