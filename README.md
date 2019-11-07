HaMStR-OneSeq - How to get started

Installation
############

HaMStR-OneSeq is a framework depending on several individual databases and data packages that can be downloaded automatically after you get your copy of the GitHub repository. Additionally there are some dependencies, which must be installed prior usage (please see Dependencies).

	1. git clone https://github.com/BIONF/hamstr

This provides you with a clone of the latest version of HaMStR. The repository does not contain the following executables (FASTA36 alignment tool set), databases, and data packages:
* fasta-36.3.8e (http://faculty.virginia.edu/wrpearson/fasta/fasta36/)
* NCBI Taxonomy, Pfam DB, SMART DB, CAST, Coils, SEG, SignalP, TMHMM.
* blast_dir, genome_dir, weight_dir

Please install the provided FASTA package (v.fasta-36.3.8e, for more information I would kindly refer you to the user manual of the FASTA package of programs).
Change into the source directory of the fasta36 package. 

	2. cd hamstr/bin/aligner/fasta-36.3.8e/src

Please select the appropriate Makefile for your system (Makefile.linux64_sse2 or Makefile.os_x86_64 are the most common cases).
With the following make command you will create the executable programs in ../bin/fasta36.
For Linux users:

	3. make -f ../make/Makefile.linux64_sse2 all

OR for MacOSX users:

	3. make -f ../make/Makefile.os_x86_64 all

Please try the following command for a quick test of your installation:

	4. ../bin/fasta36 -q ../seq/mgstm1.aa ../seq/prot_test.lseg

You can download and extract the provided data package (data_HaMStR.tar version2018 (1.7GB)) with the help of the script install_data.sh.
Please change into your main hamstr directory and follow the stepwise instructions below.

	5. cd to/your/hamstr
	6. bash ./bin/install_data.sh

Please configure your HaMStR to your personal needs and preferences with the following steps.
To set the ENV variable $ONESEQDIR add the following line to your .bashrc or .bash_profile (and reload your terminal).

	7. export ONESEQDIR=/absolute/path/to/your/hamstr

Run the provided configure script to wrap everything up. Please change into your hamstr/bin directory.
For Linux users:

	8. bash ./configure -p -n

OR for MaxOSX users:

	8. bash ./configure_mac -p -n

Your personal copy of HaMStR should be good to go.


Usage
#####

HaMStR will run smoothly with the provided sample file if everything is set correctly:

* hamstr/data/infile.fa (your input files should be placed here)
* running oneSeq.pl for the first time may take a while due to indexing steps of the NCBI taxonomy files. 

	1. `oneSeq.pl -h` (gives you an overview about available options)
	2. `oneSeq.pl -sequence_file=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global`

HaMStR-OneSeq integrates the prediction of orthologs and the comparison of their Feature Architecture Similarity score (FAS). The provided output file seqname.extended.profile contains the ID’s and the FAS scores for the orthologs in a TAB-separated file format (FAS scores are computed pairwise between the seeding input gene and it’s predicted orthologous genes). 


Visualization
#############

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


Gene sets, Annotations, Blast DBs
#################################

Within the data package (https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) we provide a set of 78 reference taxa (gene sets in genome_dir, annotations in weight_dir, blast databases in blast_dir). They can be automatically downloaded und put into place with the install_data.sh script (please see above: Installation 5. and 6.). This data comes "ready to use" with the HaMStR-OneSeq framework. Species data must be present in the three directories listed below. For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version]).:

* I.	genome_dir (Contains sub-directories for proteome fasta files for each species)
* II.	blast_dir (Contains sub-directories for BLAST databases made with makeblastdb out of your proteomes)
* III.	weight_dir (Contains sub-directories for feature annotation files for each proteome)


However, if needed the user can manually add further gene sets (multifasta format) and place them into the respective directories (genome_dir, weight_dir, blast_dir). Please note, that every taxon/species must be present in the NCBI taxonomy. The following steps need to be conducted:

* 1.) Download the gene set of your taxon of interest as amino acid sequences from the NCBI database.
* 2.) Rename the file in accordance to the naming schema of hamstr:     SPECIES@12345@1.fa
     - ([Species acronym]@[NCBI ID]@[Proteome version])
* 3.) Fasta header must be whitespace free and unique within the gene set (short header make your life easier for downstream analysis).
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

4.) After your gene set (proteomic data) is prepared and placed into the respective sub-directory in the genome_dir directory you can conduct the following instructions:

5.1) Create a Blast DB for the species within the blast_dir

	makeblastdb -dbtype prot -in genome_dir/SPECI@00001@1/SPECI@00001@1.fa -out blast_dir/SPECI@00001@1/SPECI@00001@1
	
5.2) Create a symbolic link with the blast_dir (change into the respective sub-directory in the blast_dir)
	
	cd blast_dir/SPECI@00001@1 
	ln -s ../../genome_dir/SPECI@00001@1/SPECI@00001@1.fa SPECI@00001@1.fa
	
6.) Create the annotation files for your taxon with the provided perl script

	perl /path/to/your/hamstr/bin/fas/annotation.pl -fasta=/path/to/your/hamstr/genome_dir/SPECI@00001@1/SPECI@00001@1.fa -path=/path/to/your/hamstr/weight_dir -name=SPECI@00001@1

Please take care that all parameter paths are provided as absolute paths. This action takes considerably longer than the BLAST database creation with makeblastdb (it takes about one hour to annotate a gene set with 5000 sequences).

To prove if your manually added species is integrated into the HaMStR framework your can run:

	perl bin/oneSeq.pl -showTaxa
	
This command simply prints a list of all available taxa.


Happy HaMStRing.

For further support or bug reports please contact: bergmann@bio.uni-frankfurt.de


Dependencies
* blastp (blastall)
* genewise
* hmmsearch
* hmmbuild
* clustalw2
* mafft-linsi
* Perl
* Python (FAS scoring compatible with both python v2.7 and 3.6)
    
    locale-gen
    en_US
    en_US.UTF-8
    de_DE
    de_DE.UTF-8
    
    Packages for Ubuntu 16.04:
    blast2
    hmmer
    clustalw
    mafft
    libdbi-perl
    libipc-run-perl
    wise
    locales
    ncbi-blast+
    ncbi-blast+-legacy
    lib32ncurses5
    lib32z1

    Packages for Perl:
    Parallel/ForkManager.pm
    IPC/Run.pm




