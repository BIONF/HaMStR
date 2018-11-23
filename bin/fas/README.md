# FAS - Feature Architecture Similarity
FAS is a new release of the original [FACT](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-417) algorithm. It calculates the so called FAS-score which is a measure of how similar the feature architectures of two proteins are. This is done by combining the Multiplicity Score (MS) and the Positional Score (PS) from FACT. Unlike the original FACT, FAS can resolve feature architectures that have overlapping features by searching for the best overlap-free path. This can be done either extensively or by using the priority mode, a greedy approach. FAS also allows for more options in the weighting of features. 

The main FAS script is written in Python and should run on both Python 2 and Python 3. The additional annotation script that generates the standart input for FAS is written in Perl

# Table of Contents
* [Installation](#installation)
* [Usage](#usage)
  * [How to get started](##how-to-get-started)
  * [Weighting](##weighting)
  * [In/Output Options](##in/output-options)
  * [Threshold Options](##threshold-options)
  * [Additional Options](##additional-options)
* [Contact](#contact)



# Installation
FAS is part of the [HaMStR oneseq](https://github.com/BIONF/HaMStR) package. To get it simply clone the HaMStR repository and follow the installation instructions.

```
git clone https://github.com/BIONF/HaMStR
```

Afterwards you might want to check if the Pfam HMMs under HaMStR/bin/fas/Pfam/Pfam-hmms/ are up-to-date. If not, download the newest [release](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases) and replace the corresponding files in HaMStR/bin/fas/Pfam/Pfam-hmms/.


# Usage
## How to get started
FAS comes with two main scripts: annotation.pl, which generates the standart input for FAS, and the actual FAS script greedyFAS.
To get started using FAS you need the protein sequence of the two (or more) proteins you want to compare. You should have two file in (Multi-)Fasta format, one for the seed protein(s) and one for the ortholog(s). Begin by using the annotation.pl script in the FAS directory to annotate features onto the protein sequences: 

```
perl annotation.pl -fasta=INPUTPATH/example.fa -path=OUTPUTPATH/ -name=example
```

This should give you an output folder of the chosen name containing seven xml files, one for each feature type used in the default set from FACT. Once you have annotated the features of both, the seed and ortholog proteins, you are ready to use the actual FAS algorithm with the two output folders of annotation script:

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME 
```

This will give two xml files. The first file (JOBNAME) contains the scores and chosen paths for each protein pairing between the seed and ortholog set: 

```
<?xml version="1.0"?>
<out FAS_version="1.6.0" weighting="uniform" constraints="none" MS="new" priority="50/5000" overlap="0/0.4" efilters="0.001/0.01" timelimit="off" scoreweights="(0.7, 0.0, 0.3)" cores="1" linearized="pfam/smart" normal="cast/coils/seg/signalp/tmhmm">
    <query id="EXAMPLE" length="439">
        <template id="EXAMPLE" score="1.0" MS="1.0" PS="1.0" CS="0.0" LS="1.0" length="439" mode="exhaustive/exhaustive">
            <template_path>
                <feature type="Pfam_PF00349">
                    <instance start="98" end="238"/>
                </feature>
            </template_path>
            <query_path>
                <feature type="Pfam_PF00349">
                    <instance start="98" end="238"/>
                </feature>
            </query_path>
        </template>
    </query>
</out>
```

The root node (out) contains information about the run: First the FAS version used, then the weighting that was used (if a reference was given the correction function will be added, for example: /loge), the name of the constraint file if given, the MS (multiplicity score) scheme used, the two priority thresholds (feature number/cardinality), the overlap thresholds (size/percentage), the two e-filters (type/instance), the timelimit, the weights of the three scores (MS, CS, PS), the number of cores used, the linearized features and finally the normal features.

The root node has one type of child node (query) which is the proteins from the ortholog set (one node for each protein). It contains the id of the protein and the sequence length. Each query node has one type of child node, the template which is the proteins from the seed set. This node also contains the id and length of the protein together with all scores and the mode with which the path was chosen for the seed (1.) and query (2.). There are two possibilities for mode: exhaustive will always give the best path, priority might give a non optimal solution. The template node always has exactly two child nodes, the template path and the query path, which contain the chosen paths for both proteins for the score calculation. They have one type of child node (feature) which contains the feature id/name (not that the inputfile name is added at the front, like: Pfam_) and, if a reference was given, its weight in the score. If no reference was given, all features are weighted equally. Each feature node has a number of instance child nodes which list all occurrences of the feature in the path with their start and stop position. 

The second file (JOBNAME_architecture) contains the architecture of each protein from both sets: 

```
<?xml version="1.0"?>
<architectures>
    <template id="EXAMPLE" length="439">
        <architecture>
            <feature type="Pfam_PF00349" evalue="NULL">
                <instance inst_eval="NULL" start="98" end="238"/>
            </feature>
        </architecture>
    </template>
    <query id="EXAMPLE" length="439">
        <architecture>
            <feature type="Pfam_PF00349" evalue="NULL">
                <instance inst_eval="NULL" start="98" end="238"/>
            </feature>
        </architecture>
    </query>
</architectures>
```

Above we used FAS at its very basic with no more options given then the bare necessities. However, there are a lot of options that can be changed. Remeber that you can always use the help option (-h) to display a description of all options. 

## Weighting
In the first example we used a uniform weighting, so that all features would be weighted equally. To use a different weighting we need give a reference protein set (-r): 

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -e PATH/reference 
```

The reference set is there give FAS information on the abundance of the features so that it can create a weighting based on that. It should have the same format as the two input sets (seed&ortholog). If you are not sure what set to use you should go with the proteome of the species the ortholog protein came from. 

You can use the weight correction option (-g) to change the correction function on the feature counts (default: loge): 

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -r PATH/reference -g linear 
```

This sets the weighting to the original from FACT.

Finally you can give FAS a file with weight constraints (-x): 

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -r PATH/reference -x PATH/weight_constraints 
```

This file should look like this: 

```
 #tool_constraints (this line must be kept)
 cast N
 coils N
 pfam 0.5
 seg 0.1
 signalp N
 smart N
 tmhmm N
 #feature_constraints (this line must be kept)
 tmhmm_transmembrane 0.1
 pfam_HlyD 0.25
```

The sum of all constraint values under #tool_constraints as well as #feature_constraints should not exceed 1.0. 
Finally you can also change the weighting of the individual scores themselves with the weights option (-w): 

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -w (0.6, 0.1, 0.3) 
```

This would set the weight of the MS to 0.6, the CS to 0.1 and the PS to 0.3. The CS will only be calculated if its weight is higher than 0. (Not recommended to change this)

## In/Output Options
You can set the extendedout option (-e) to 0 to deactivate the _architecture output: 

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -e 0
```

You can use the raw_output option (-a) to make FAS output its score to STDOUT by setting it to 1 or 2. 1 also suppresses the normal xml output: 

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -a 1 
```

Using the feature_info (-y) which gives you a file with information on the abundance of all seed and query features in the reference

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -y 1 
```

Finally you can use the featuretypes option (-d) to change what feature file FAS will use and which will be linearized:

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -d PATH/featuretypes
```

This needs a file that looks like this:

```
 #linearized
 pfam
 smart
 #normal
 cast
 coils
 seg
 signalp
 tmhmm
 newfeatures
```

The features under #linearized will get linearized while the ones under #normal are allowed to overlap. Note, that FAS will search for files with the names given here (pfam.xml, smart.xml, cast.xml,…, newfeatures.xml). All these files need to be in one of the xml input formats of FAS: 
```
<?xml version="1.0"?>
<tool name="newfeatures">
 <protein id="EXAMPLE" length="439">
     <feature type="cd00012" instance="1">
         <start start="106"/>
         <end end="270"/>
     </feature>
 </protein>
</tool>
```

OR

```
<?xml version="1.0"?>
<tool name="PfamDB">
    <protein id="EXAMPLE" length="439">
        <feature type="Guanylate_cyc" instance="1" clan="CL0276" evalue="2.8e-10">
            <start inst_eval="1.8e-14" start="1" end="48"/>
        </feature>
    </protein>
</tool>
```

## Threshold Options

You can change the e-value filter for features with the efilter option (-f) and for feature instances with the inst_efilter (-i):

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -f 0.002 -i 0.02
```

You can also change the thresholds for priority mode with the priority_threshold option (-t) for the total number of feature instances and the max_cardinality option (-m) for the number of paths:

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -t 40 -m 4000 
```

There is also the possibility of allowing small overlaps in the path. This can be done by setting max_overlap (-c) to an overlap size (in amino acids) of your choice:

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -c 10
```

By default overlaps are allowed to take up at maximum 40% of the size of either feature. This can be changed by using the max_overlap_percentage (set by using a value between 0.0 and 1.0).

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME -c 10 --max_overlap_percentage 0.5
```

This would set the threshold to 50%.

## Additional Options

The prior checks using the priority thresholds can be deactivated using the priority_check option (not recommended):

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME --priority_check 0
```

This opens up two more options: First the timelimit option. This option allows you to set(in seconds)/deactivate(set to 0) the timelimit the calculation for each pair is allowed to take before priority mode takes over (default 7200s). Secondly, the cores option which allows you to run the calculation on multiple cores:

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME --priority_check 0 --cores 4 --timelimit 600
```

You can switch between classic and new MS calculation using the classicMS option:

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME --classicMS 0
```

Finally, the bidirectional option (-b) tells FAS to run the scoring in both directions if you set it to 1, you can also set this to 2 to use to cores for this(should not be used together with the cores option).

```
python greedyFAS.py -q PATH/ortholog -s PATH/seed -j JOBNAME --bidirectional 1
``` 

This will generate two normal output-files. One for each scoring direction. Additionally, it creates a csv file that gives a short overview over all scores.



# Contact
dosch@bio.uni-frankfurt.de
