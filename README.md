# HaMStR-OneSeq
[![PyPI version](https://badge.fury.io/py/h1s.svg)](https://pypi.org/project/h1s/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/BIONF/HaMStR.svg?branch=master)](https://travis-ci.com/BIONF/HaMStR)

# Table of Contents
* [How to install](#how-to-install)
     * [Install the h1s package](#install-the-h1s-package)
     * [Setup HaMStR-oneSeq](#setup-hamstr-oneseq)
* [Usage](#usage)
* [HaMStR-oneSeq data set](#hamstr-oneseq-data-set)
     * [Adding a new gene set into HaMStR-oneSeq](#adding-a-new-gene-set-into-hamstr-oneseq)
     * [Adding a list of gene sets into HaMStR-oneSeq](#adding-a-list-of-gene-sets-into-hamstr-oneseq)
* [Bugs](#bugs)
* [How to cite](#how-to-cite)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*HaMStR-oneSeq* is distributed as a python package called *h1s*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

## Install the h1s package
You can install *h1s* using `pip`:
```
python3 -m pip install h1s
```

or, in case you do not have admin rights, and don't use package systems like Anaconda to manage environments you need to use the `--user` option:
```
python3 -m pip install --user h1s
```

and then add the following line to the end of your **~/.bashrc** or **~/.bash_profile** file, restart the current terminal to apply the change (or type `source ~/.bashrc`):

```
export PATH=$HOME/.local/bin:$PATH
```

## Setup HaMStR-oneSeq

After installing *h1s*, you need to setup *HaMStR-oneSeq* to get its dependencies and pre-calculated data.

You can do it by just running this command
```
setup1s -o /output/path/for/oneSeq/data
```
or, in case you are using Anaconda
```
setup1s -o /output/path/for/oneSeq/data --conda
```

*You should have the sudo password ready, otherwise some missing dependencies cannot be installed. See [dependency list](#dependencies) for more info. If you do not have root privileges, ask your admin to install those dependencies using `setup1s --lib` command.*

[Pre-calculated data set](https://github.com/BIONF/HaMStR/wiki/Input-and-Output-Files#data-structure) of HaMStR-oneSeq will be saved in `/output/path/for/oneSeq/data`. After the setup run successfully, you can start using *HaMStR-oneSeq*.

*For debugging the setup, please create a log file by running the setup as e.g. `setup1s | tee log.txt` for Linux/MacOS or `setup1s --conda | tee log.txt` for Anaconda and send us that log file, so that we can trouble shoot the issues. Most of the problems can be solved by just re-running the setup.*

# Usage
*HaMStR-oneSeq* will run smoothly with the provided sample input file in 'infile.fa' if everything is set correctly.

```
h1s --seqFile infile.fa --seqName test --refspec HUMAN@9606@3
```
The output files with the prefix `test` will be saved at your current working directory.
You can have an overview about all available options with the command
```
h1s -h
```

Please find more information in [our wiki](https://github.com/BIONF/HaMStR/wiki) to learn about the [input and outputs files](https://github.com/BIONF/HaMStR/wiki/Input-and-Output-Files) of *HaMStR-oneSeq*.

# HaMStR-oneSeq data set

Within the data package we provide a set of 78 reference taxa. They can be automatically downloaded during the setup. This data comes "ready to use" with the *HaMStR-OneSeq* framework. Species data must be present in the three directories listed below:

* genome_dir (Contains sub-directories for proteome fasta files for each species)
* blast_dir (Contains sub-directories for BLAST databases made with `makeblastdb` out of your proteomes)
* weight_dir (Contains feature annotation files for each proteome)

For each species/taxon there is a sub-directory named in accordance to the naming schema ([Species acronym]@[NCBI ID]@[Proteome version])

HaMStR-oneSeq is not limited to those 78 taxa. If needed the user can manually add further gene sets (multifasta format) using provided python scripts.

## Adding a new gene set into HaMStR-oneSeq
For adding **one gene set**, please use the `addTaxon1s` function:
```
addTaxon1s -f newTaxon.fa -i tax_id [-o /output/directory] [-n abbr_tax_name] [-c] [-v protein_version] [-a]
```

in which, the first 3 arguments are required including `newTaxon.fa` is the gene set that need to be added, `tax_id` is its NCBI taxonomy ID, `/output/directory` is where the sub-directories can be found (*genome_dir*, *blast_dir* and *weight_dir*). If not given, new taxon will be added into the same directory of pre-calculated data. Other arguments are optional, which are `-n` for specify your own taxon name (if not given, an abbriviate name will be suggested based on the NCBI taxon name of the input `tax_id`), `-c` for calculating the BLAST DB (only needed if you need to include your new taxon into the list of taxa for compilating the core set), `-v` for identifying the genome/proteome version (default will be 1), and `-a` for turning off the annotation step (*not recommended*).

## Adding a list of gene sets into HaMStR-oneSeq
For adding **more than one gene set**, please use the `addTaxa1s` script:
```
addTaxa1s -i /path/to/newtaxa/fasta -m mapping_file [-o /output/directory] [-c]
```
in which, `/path/to/taxa/fasta` is a folder where the FASTA files of all new taxa can be found. `mapping_file` is a tab-delimited text file, where you provide the taxonomy IDs that stick with the FASTA files:

```
#filename	tax_id	abbr_tax_name	version
filename1.fa	12345678
filename2.faa	9606
filename3.fasta	4932	my_fungi
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
