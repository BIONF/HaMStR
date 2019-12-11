# HaMStR-OneSeq

## Installation

### 0. Basic system tools requirement
You need to have `wget`, `grep` and `sed` (or `gsed` for **MacOS**) to install HaMStR. So please install them if they are missing. For MacOS users, we recommend using [Homebrew](https://brew.sh) to install those command line tools.

### 1a. Install in Ubuntu/MacOS

Get HaMStR source code from GitHub
```
git clone --depth=1 https://github.com/BIONF/HaMStR
```

Run `setup.sh` script in the HaMStR/bin folder to install HaMStR and its dependencies
```
cd HaMStR
bin/setup.sh
```
*Enter root password if required (some dependencies need root privileges to be installed. See [dependency list](#dependencies) for more info.)*

After the setup run successfully, you can start using HaMStR (in some cases you should restart the terminal).

### 1b. Install using Anaconda

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
conda install -c trvinh hamstr
setup_hamstr
```

After the setup run successfully, you can start using HaMStR (in some cases you should restart the terminal).

*For debugging the installation, please create a log file by running the setup as e.g. `bin/setup.sh | tee log.txt` for Linux/MacOS or `setup_hamstr | tee log.txt` for Anaconda and send us that log file, so that we can trouble shoot the issues. Most of the problems can be solved by just re-running the setup.*

## Usage
HaMStR will run smoothly with the provided sample file if everything is set correctly:

* hamstr/data/infile.fa (your input files should be placed here)
* running oneSeq.pl for the first time may take a while due to indexing steps of the NCBI taxonomy files.

	1. `oneSeq.pl -h` (gives you an overview about available options)
	2. `oneSeq.pl -sequence_file=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global`

HaMStR-OneSeq integrates the prediction of orthologs and the comparison of their Feature Architecture Similarity score (FAS). The provided output file seqname.extended.profile contains the ID’s and the FAS scores for the orthologs in a TAB-separated file format (FAS scores are computed pairwise between the seeding input gene and it’s predicted orthologous genes).


## Visualization
HaMStR provides you with a set of output files in plain text format. For a rich visualisation of the provided information you can plug them into the Phyloprofile tool (https://github.com/BIONF/phyloprofile)

* seqname.extended.profile  (Profile with FAS score for ortholgs)
* seqname.extended.fa  (Orthologs file in Fasta format)
* seqname_1.matrix  (Input file for Phyloprofile tool)
* seqname_1.domains  (Additional input file for Phyloprofile tool)
* seqname_0.matrix  (Input file for Phyloprofile, if opt -countercheck is set)
* seqname_0.matrix  (Additional input file for Phyloprofile, if opt -countercheck is set)

You can combine multiple HaMStR-OneSeq runs into a single Matrix/phylogenetic profile for data visualisation and data exploration. Each run is identified by the given seqname (opt -seqname=<>). This is either given by the user or randomly assigned. The following steps are necessary:

Concatenate all desired profile files (*.profile) into one combined profile:

	1. cat seqname1.extended.profile seqname2.extended.profile seqname3.extended.profile > combined.profile

Re-run the parsing script provided in HaMStR/bin/visuals/parseOneSeq.pl from your current data directory with the combined profile.

	2. perl /path/to/HaMStR/bin/visuals/parseOneSeq.pl -i combined.profile -o combined

This provides you with a combined matrix file (combined_1.matrix), which is suitable as input for Phyloprofile. If you have a two-way profile containing the forward (1) and the backward (0) FAS scores you get a combined matrix for both directions.

* combined_1.matrix
* combined_0.matrix (optionally, trigged in your oneSeq.pl command with the option -countercheck)

To prepare the additional input file (*.domains) you just need to concatenate them with each other (please mind the distinction between forward (1) and backward (0) comparisons and do not mix them up).

	3. cat seqname1_1.domains seqname2_1.domains seqname3_1.domains > combined_1.domains

The resulting file combined_1.matrix and combined_1.domains can be plugged into the Phyloprofile tool (R shiny) for further investigation.


## Gene sets, Annotations, Blast DBs

Within the data package (https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) we provide a set of 78 reference taxa (gene sets in genome_dir, annotations in weight_dir, blast databases in blast_dir). They can be automatically downloaded with the `setup.sh` script. This data comes "ready to use" with the HaMStR-OneSeq framework. Species data must be present in the three directories listed below. For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version]).:

* genome_dir (Contains sub-directories for proteome fasta files for each species)
* blast_dir (Contains sub-directories for BLAST databases made with makeblastdb out of your proteomes)
* weight_dir (Contains sub-directories for feature annotation files for each proteome)


However, if needed the user can manually add further gene sets (multifasta format) and place them into the respective directories (genome_dir, weight_dir, blast_dir). Please note, that every taxon/species must be present in the NCBI taxonomy. The following steps need to be conducted:

1) Download the gene set of your taxon of interest as amino acid sequences from the NCBI database.

2) Rename the file in accordance to the naming schema of hamstr:     SPECIES@12345@1.fa ([Species acronym]@[NCBI ID]@[Proteome version])

3) Fasta header must be whitespace free and unique within the gene set (short header make your life easier for downstream analysis).
     - the following bash command uses sed to cut the header at the first whitespace: sed -i "s/ .*//" SPECIES@12345@1.fa
     - example:

before:

	>EXR66326.1 biofilm-associated domain protein, partial [Acinetobacter baumannii 339786]
	MTGEGPVAIHAEAVDAQGNVDVADADVTLTIDTTPQDLITAITVPEDLNGDGILNAAELGTDGSFNAQVALGPDAVDGTV
	VNVNGTNYTVTAADLANGYITATLDATAADPVTGQIVIHAEAVDAQGNVD
	>EXR66351.1 hypothetical protein J700_4015, partial [Acinetobacter baumannii 339786]
	NRRLLITTQPTATDSNYKTPIYINAPNGELYFANQDETSVSSVVFKRVIGATAANAPYVASDSWTKKIRKWNTYNHEVSK
	VGRFIAPMMLTYDVTFTTQQNNAGWSISKESTGVYRLQRDSGVTTELANPHIEVSGIFAGTGLGSGDVILPPTLQAIEAY
	>EXR66376.1 bacterial Ig-like domain family protein, partial [Acinetobacter baumannii 339786]
	DGVDYPAVNNGDGTWTLADNTLPTLADGPHTITVTATDAAGNVGNDTAVVTIDTVAPNAPVLDPINATDPVSGQAEPGST
	VTVTYPDGTTATVVAGTDGSWSVPNPGNLVDGDTVTATAT
	...
after (this is how your sequence data should look like):

	>EXR66326.1
	MTGEGPVAIHAEAVDAQGNVDVADADVTLTIDTTPQDLITAITVPEDLNGDGILNAAELGTDGSFNAQVALGPDAVDGTV
	VNVNGTNYTVTAADLANGYITATLDATAADPVTGQIVIHAEAVDAQGNVD
	>EXR66351.1
	NRRLLITTQPTATDSNYKTPIYINAPNGELYFANQDETSVSSVVFKRVIGATAANAPYVASDSWTKKIRKWNTYNHEVSK
	VGRFIAPMMLTYDVTFTTQQNNAGWSISKESTGVYRLQRDSGVTTELANPHIEVSGIFAGTGLGSGDVILPPTLQAIEAY
	>EXR66376.1
	DGVDYPAVNNGDGTWTLADNTLPTLADGPHTITVTATDAAGNVGNDTAVVTIDTVAPNAPVLDPINATDPVSGQAEPGST
	VTVTYPDGTTATVVAGTDGSWSVPNPGNLVDGDTVTATAT

4) After your gene set (proteomic data) is prepared and placed into the respective sub-directory in the genome_dir directory you can conduct the following instructions:

5.1) Create a Blast DB for the species within the blast_dir

	makeblastdb -dbtype prot -in genome_dir/SPECI@00001@1/SPECI@00001@1.fa -out blast_dir/SPECI@00001@1/SPECI@00001@1

5.2) Create a symbolic link with the blast_dir (change into the respective sub-directory in the blast_dir)

	cd blast_dir/SPECI@00001@1
	ln -s ../../genome_dir/SPECI@00001@1/SPECI@00001@1.fa SPECI@00001@1.fa

6) Create the annotation files for your taxon with the provided perl script

	perl /path/to/your/hamstr/bin/fas/annotation.pl -fasta=/path/to/your/hamstr/genome_dir/SPECI@00001@1/SPECI@00001@1.fa -path=/path/to/your/hamstr/weight_dir -name=SPECI@00001@1

Please take care that all parameter paths are provided as absolute paths. This action takes considerably longer than the BLAST database creation with makeblastdb (it takes about one hour to annotate a gene set with 5000 sequences).

To prove if your manually added species is integrated into the HaMStR framework your can run:

	perl bin/oneSeq.pl -showTaxa
This command simply prints a list of all available taxa.

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

_**Note: After having all these dependencies installed, you still need to run the setup script to configure HaMStR!!!**_

### System tools/libraries
* grep (ggrep)
* sed (gsed)
* wget (wget)
* build-essential
* curl (curl)
* locales
* lib32ncurses5
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

## Cite
Ebersberger, I., Strauss, S. & von Haeseler, A. HaMStR: Profile hidden markov model based search for orthologs in ESTs. BMC Evol Biol 9, 157 (2009) doi:10.1186/1471-2148-9-157

## Contact
For further support or bug reports please contact: ebersberger@bio.uni-frankfurt.de
