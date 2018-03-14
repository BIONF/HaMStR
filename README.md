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

You can download and extract the provided data package (data_HaMStR.tar) with the help of the script install_data.sh.
Please change into your main hamstr directory and follow the stepwise instructions below.

	5. cd to/your/hamstr
	6. bash ./bin/install_data.sh

Please configure your HaMStR to your personal needs and preferences with the following steps.
To set the ENV variable $ONESEQDIR add the following line to your .bashrc or .bash_profile.

	7. export ONESEQDIR=/absolute/path/to/your/hamstr

Run the provided configure script to wrap everything up.
For Linux users:

	8. bash ./bin/configure -p -n

OR for MaxOSX users:

	8. bash ./bin/configure_mac -p -n

Your personal copy of HaMStR should be good to go.

Usage
#####

HaMStR will run smoothly with the provided sample file if everything is set correctly:

* hamstr/data/infile.fa (your input files should be placed here)

	1. oneSeq.pl -help (gives you an overview about available options)
	2. oneSeq.pl -sequence_file=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global

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
* Python (FAS scoring requires python v2.7)
    
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




