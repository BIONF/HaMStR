# HaMStR-OneSeq
[![conda-install](https://anaconda.org/bionf/hamstr/badges/installer/conda.svg)](https://anaconda.org/bionf/hamstr)
[![conda-version](https://anaconda.org/bionf/hamstr/badges/version.svg)](https://anaconda.org/bionf/hamstr)
[![GPLv3-license](https://anaconda.org/bionf/hamstr/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.de.html)

# Table of Contents
* [How to install](#how-to-install)
     * [0. Basic system tools requirement](#0-basic-system-tools-requirement)
     * [1a. Install using Anaconda](#1a-install-using-anaconda)
     * [1b. Install in Ubuntu/MacOS](#1b-install-in-ubuntumacos)
* [Usage](#usage)
* [HaMStR and the utilisation of FAS](#hamstr-and-the-utilisation-of-fas)
* [Output visualization using PhyloProfile](#output-visualization-using-phyloprofile)
* [Pre-calculated data set](#pre-calculated-data-set)
* [Dependencies](#dependencies)
  * [System tools/libraries](#system-toolslibraries)
  * [Bioinformatics tools](#bioinformatics-tools)
  * [Perl modules](#perl-modules)
* [How to cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

## How to install

### 0. Basic system tools requirement
You need to have `wget`, `grep` and `sed` (or `gsed` for **MacOS**) to install HaMStR. So please install them if they are missing. For MacOS users, we recommend using [Homebrew](https://brew.sh) to install those command line tools.
To use [FAS tool](https://github.com/BIONF/FAS) (a dependency of HaMStR), you also need [Python 3](https://www.python.org/downloads/).

### 1a. Install using Anaconda

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

### 1b. Install in Ubuntu/MacOS

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

## Usage
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

The output consist of these text files (*note: `test` is your defined -seqName parameter*)
1) `test.extended.fa`: a multiple FASTA file containing ortholog sequences and the query gene
2) `test.extended.profile`: a tab-delimited file containing list of orthologous sequences and their correspoding similarity scores by comparing their feature architectures with the one of the query gene (for more info about this score, please read [this document](https://bionf.github.io/FAS))
3) `test.phyloprofile`: an input file for visualisation the phylogenetic profile of the query gene using [PhyloProfile tool](https://github.com/BIONF/phyloprofile)
4) `test_1.domains` (and optional, `test_0.domains`): a protein domain annotation file for all the sequences present in the orthologous group. The `_0` or `_1` suffix indicates the direction of the feature architecture comparison, in which `_1` (forward) means that the query gene is used as *seed* and it orthologs as *target* for the comparison, while `_0` (backward) is vice versa.

## HaMStR and the utilisation of FAS
HaMStR integrates the prediction of orthologs and the calculation of the Feature Architecture Similarty (FAS) scores. FAS scores are computed pairwise between the query gene and it's predicted orthologous genes using [FAS tool](https://github.com/BIONF/FAS), which will be automatically installed during the setup of HaMStR.

## Output visualization using PhyloProfile
For a rich visualisation of the provided information from the HaMStR outputs, you can plug them into the [Phyloprofile tool](https://github.com/BIONF/phyloprofile).

The main input file for *PhyloProfile* is `seqname.phyloprofile`, which contains list of all orthologous gene names and the taxonomy IDs of their taxa together with the FAS scores (if available). For analysing more information such as the FASTA sequences or the domain annotations, you can optionally input `seqname.extended.fa` and `seqname_1.domains` (or `seqname_0.domains`) to *PhyloProfile*.

You can combine multiple HaMStR runs into a single phylogenetic profile input for data visualisation and data exploration. Each run is identified by the given seqname (opt -seqname=<>). This is either given by the user or randomly assigned. The following steps are necessary:

```
# concatenate all desired profile files into one combined profile

cat *.extended.profile > combined.extended.profile

# re-run the parsing script from your current data directory with the combined profile

perl /path/to/HaMStR/bin/visuals/parseOneSeq.pl -i combined.extended.profile -o combined.phyloprofile
```

To prepare the additional input file (*.domains) you just need to concatenate them with each other (please mind the distinction between forward (1) and backward (0) FAS comparisons and do not mix them up).

```
cat *_0.domains > combined_0.domains
cat *_1.domains > combined_1.domains
```

The resulting file `combined.phyloprofile`, `combined_0.matrix` and `combined_1.domains` can be then plugged into the *Phyloprofile tool* for further investigation.


## Pre-calculated data set

Within the data package (https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) we provide a set of 78 reference taxa (gene sets in genome_dir, annotations in weight_dir, blast databases in blast_dir). They can be automatically downloaded during the setup. This data comes "ready to use" with the HaMStR-OneSeq framework. Species data must be present in the three directories listed below. For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version]):

* genome_dir (Contains sub-directories for proteome fasta files for each species)
* blast_dir (Contains sub-directories for BLAST databases made with makeblastdb out of your proteomes)
* weight_dir (Contains sub-directories for feature annotation files for each proteome)

However, if needed the user can manually add further gene sets (multifasta format) using provided python scripts.

### Adding a new gene set into HaMStR
For adding **one gene set**, please use the `bin/addTaxonHamstr.py` script:
```
python3 bin/addTaxonHamstr.py -f newTaxon.fa -n abbr_name -i tax_id -o /path/to/HaMStR [-c] [-v protein_version] [-a]
```

in which, the first 3 arguments are required including `abbr_name` is the species acronym name, `tax_id` is its NCBI taxonomy ID, `/path/to/HaMStR` is where the sub-directories will be saved (genome_dir, blast_dir and weight_dir). Other arguments are optional, which are `-c` for calculating the BLAST DB (only needed if you need to include your new taxon into the list of taxa for compilating the core set), `-v` for identifying the genome/proteome version (default will be 1), and `-a` for turning off the annotation step (*not recommended*).

### Adding a list of gene sets into HaMStR
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

Please check this [wiki page](https://github.com/BIONF/HaMStR/wiki/Add-new-taxa-to-HaMStR) for more details.

## Dependencies
HaMStR has some dependencies, that either will be automatically installed via the setup script, or must be installed by your system admin if you don't have the root privileges. In the following you will find the full list of HaMStR's dependencies for Ubuntu system as well as the alternatives for MacOS. In Ubuntu, you can install those system and bioinformatics tools/libraries using `apt-get` tool
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

_**Note: After having all these dependencies installed, you still need to run the setup script to configure HaMStR!!!**_

### System tools/libraries
* grep (ggrep)
* sed (gsed)
* wget (wget)
* build-essential
* curl (curl)
* locales
* lib32ncurses6
* lib32z1

*(In parentheses are Mac's alternative tools)*

### Bioinformatics tools
* wise (brewsci/bio/genewise)
* hmmer (hmmer)
* ncbi-blast+ (blast)
* blast2
* clustalw (brewsci/bio/clustal-w)
* mafft (mafft)
* muscle (brewsci/bio/muscle)

*(In parentheses are Mac's alternative tools)*

### Perl modules
* libdbi-perl
* libipc-run-perl
* perl-doc
* DBI
* DB_File
* File::Copy
* File::Path
* File::Basename
* File::Which
* List::Util
* Parallel::ForkManager
* POSIX
* XML::SAX
* XML::NamespaceSupport
* XML::Parser
* Getopt::Long
* IO::Handle
* IPC::Run
* Statistics::R
* Term::Cap
* Time::HiRes
* Bio::AlignIO
* Bio::Align::ProteinStatistics
* Bio::DB::Taxonomy
* Bio::SearchIO
* Bio::SearchIO::blastxml
* Bio::Search::Hit::BlastHit
* Bio::Seq
* Bio::SeqIO
* Bio::SeqUtils
* Bio::Tree::Tree
* Bio::Tools::Run::StandAloneBlast

## How to cite
Ebersberger, I., Strauss, S. & von Haeseler, A. HaMStR: Profile hidden markov model based search for orthologs in ESTs. BMC Evol Biol 9, 157 (2009), [doi:10.1186/1471-2148-9-157](https://doi.org/10.1186/1471-2148-9-157)

# Contributors
- [Ingo Ebersberger](https://github.com/ebersber)
- [Vinh Tran](https://github.com/trvinh)
- [Holger Bergmann](https://github.com/holgerbgm)

## Contact
For further support or bug reports please contact: ebersberger@bio.uni-frankfurt.de
