HaMStR-OneSeq - How to get started

Installation

HaMStR-OneSeq is a framework depending on several individual databases and data packages that can be downloaded automatically after you get your copy of the GitHub repository. Additionally there are some dependencies, which must be installed prior usage (please see Dependencies).

    git clone https://github.com/BIONF/hamstr

This provides you with a clone of the latest version of HaMStR. The repository does not contain the following databases and data packages:

    NCBI Taxonomy, Pfam DB, SMART DB, CAST, Coils, SEG, SignalP, TMHMM.
    blast_dir, genome_dir, weight_dir

To download the provided data (data_HaMStR.tar):

    cd HaMStR
    bash ./bin/install_data.sh

Please configure your HaMStR to your personal needs and preferences with the following steps.
To set the ENV variable $ONESEQDIR add the following line to your .bashrc or .bash_profile

    export ONESEQDIR=/absolute/path/to/HaMStR


Run the provided configure script to wrap everything up:

    bash ./bin/configure -p -n

Your personal copy of HaMStR should good to go.

Usage

HaMStR will run smoothly with the provided sample file if everything is set correctly:

    HaMStR/data/infile.fa (your input files should be placed here)


    oneSeq.pl -help (gives you an overview about available options)
    oneSeq.pl -sequence_file=infile.fa -seqid=P83876 -refspec=HUMAN@9606@1 -minDist=genus -maxDist=kingdom -coreOrth=5 -cleanup -global

HaMStR-OneSeq integrates the prediction of orthologs and the comparison of their Feature Architecture Similarity score (FAS). The provided output file seqname.extended.profile contains the ID’s and the FAS scores for the orthologs in a TAB-separated file format (FAS scores are computed pairwise between the seeding input gene and it’s predicted orthologous genes). 


Visualization

HaMStR provides you with a set of output files in plain text format. For a rich visualisation of the provided information you can plug them into the Phyloprofile tool (https://github.com/BIONF/phyloprofile)

    seqname.extended.profile  (Profile with FAS score for ortholgs)
    seqname.extended.fa  (Orthologs file in Fasta format)
    seqname_1.matrix  (Input file for Phyloprofile tool)
    seqname_1.domains  (Additional input file for Phyloprofile tool)
    seqname_0.matrix  (Input file for Phyloprofile, if opt -countercheck is set)
    seqname_0.matrix  (Additional input file for Phyloprofile, if opt -countercheck is set)


You can combine multiple HaMStR-OneSeq runs into a single Matrix/phylogenetic profile for data visualisation and data exploration. Each run is identified by the given seqname (opt -seqname=<>). This is either given by the user or randomly assigned. The following steps are necessary:

Concatenate all desired profile files (*.profile) into one combined profile:

    cat seqname1.extended.profile seqname2.extended.profile seqname3.extended.profile > combined.profile

Re-run the parsing script provided in HaMStR/bin/visuals/parseOneSeq.pl from your current data directory with the combined profile.

    perl /path/to/HaMStR/bin/visuals/parseOneSeq.pl -i combined.profile -o combined

This provides you with a combined matrix file (combined_1.matrix), which is suitable as input for Phyloprofile. If you have a two-way profile containing the forward (1) and the backward (0) FAS scores you get a combined matrix for both directions.

    combined_1.matrix
    combined_0.matrix (optionally, trigged in your oneSeq.pl command with the option -countercheck)


To prepare the additional input file (*.domains) you just need to concatenate them with each other (please mind the distinction between forward (1) and backward (0) comparisons and do not mix them up).

    cat seqname1_1.domains seqname2_1.domains seqname3_1.domains > combined_1.domains

The resulting file combined_1.matrix and combined_1.domains can be plugged into the Phyloprofile tool (R shiny) for further investigation.


Happy HaMStRing.

For further support or bug reports please contact: bergmann[at]bio.uni-frankfurt.de


Dependencies

    blastall
    genewise
    hmmsearch
    hmmbuild
    clustalw2
    mafft-linsi
    Perl
    Python
    
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




