########################
 hamstr1seq distribution

###########
1. Overview

1.1 Directories:
	bin
	blast_dir (initially empty dir, files not included in the repository)
	core_orthologs (initially empty dir)
	data (initially empty dir)
	genome_dir (initially empty dir, files not included in the repository)
	lib
	output (initially empty dir)
	taxonomy (initially empty dir, files not included in the repository)
	tmp (initially empty dir)
	weight_dir (initially empty dir, files not included in the repository)

	Files:
	BLOSUM62.txt
	
1.2 bin directory
	configure
	configure_mac
	hamstr.pl
	multifasta2oneSeq.pl
	nentferner.pl
	oneSeq.pl
	run-query.sh
	translate.pl

1.2.1 bin/fas
	CAST (lib for: 32bit architecture support required)
	COILS2
	Pfam (not included in the repository)
	SEG
	SignalP
	SMART
	TMHMM
	
	Files:
	annotation.pl
	greedyFAS.py
1.2.2 bin/visuals
	parseArchitecture.pl	(parse FAS output files to create domain files for phyloprofile)
	parseOneSeq.pl	(parse HaMStR-OneSeq output to create matrix file for phyloprofile)
1.2.3 bin/aligner
	fasta-36.3.8e (FASTA suite of programs: https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml)
	README

1.3 Installation and Usage
	Please see README.md
