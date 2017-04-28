### README ###
Universal Protein Resource (UniProt)
====================================


The Universal Protein Resource (UniProt), a collaboration between the European
Bioinformatics Institute (EBI), the SIB Swiss Institute of Bioinformatics, and
the Protein Information Resource (PIR), is comprised of three databases, each
optimized for different uses. The UniProt Knowledgebase (UniProtKB) is the
central access point for extensively curated protein information, including
function, classification and cross-references. The UniProt Reference Clusters
(UniRef) combine closely related sequences into a single record to speed up
sequence similarity searches. The UniProt Archive (UniParc) is a comprehensive
repository of all protein sequences, consisting only of unique identifiers and
sequences.


Reference Proteomes QfO release
===============================

Some proteomes have been (manually and algorithmically) selected as reference
proteomes. They cover well-studied model organisms and other organisms of
interest for biomedical research and phylogeny.

Based on UniProt Release 2017_04, Ensembl release 87 and Ensembl Genome
release 34

This folder, ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/,
contains the following files, updated once a year:

- README (this file)

- Taxonomic division folders (Archaea, Bacteria and Eukaryota) containing
  *.gene2acc, *.fasta (canonical sequences), *_DNA.fasta, *.xml and *.idmapping
  files for 78 species, plus 42 *_additional.fasta (isoform sequences) and 44
  *_DNA.miss files

- QfO_release_2017_04.tar.gz is an archive that contains all files above.


Release 2017_04, 12-Apr-2017

Statistics (Total of species = 78):

    #(1) Number of entries in main fasta (canonical)
    #(2) Number of entries in additional fasta (isoforms)
    #(3) Number of entries in gene2acc mapping file

Proteome_ID Tax_ID  OSCODE     #(1)    #(2)    #(3)  Species Name
UP000007062 7165    ANOGA     11930    1442   14384  Anopheles gambiae (African malaria mosquito)
UP000000798 224324  AQUAE      1553       0    1557  Aquifex aeolicus (strain VF5)
UP000006548 3702    ARATH     27502   13831   41673  Arabidopsis thaliana (Mouse-ear cress)
UP000001570 224308  BACSU      4197       7    4205  Bacillus subtilis (strain 168)
UP000001414 226186  BACTN      4782       0    4823  Bacteroides thetaiotaomicron (strain ATCC 29148 / DSM 2079 / NCTC 10582 / E50 / VPI-5482)
UP000007241 684364  BATDJ      8610       0    8700  Batrachochytrium dendrobatidis (strain JAM81 / FGSC 10211) (Frog chytrid fungus)
UP000009136 9913    BOVIN     21987    2529   24612  Bos taurus (Bovine)
UP000002526 224911  BRADU      8253       0    8314  Bradyrhizobium diazoefficiens (strain JCM 10833 / IAM 13628 / NBRC 14792 / USDA 110)
UP000001554 7739    BRAFL     28542       2   28639  Branchiostoma floridae (Florida lancelet) (Amphioxus)
UP000001940 6239    CAEEL     20057    7998   28264  Caenorhabditis elegans
UP000000559 237561  CANAL      6153       1   10852  Candida albicans (strain SC5314 / ATCC MYA-2876) (Yeast)
UP000002254 9615    CANLF     20141    5395   25631  Canis lupus familiaris (Dog) (Canis familiaris)
UP000000431 272561  CHLTR       895       0     895  Chlamydia trachomatis (strain D/UW-3/Cx)
UP000006906 3055    CHLRE     14271      67   14483  Chlamydomonas reinhardtii (Chlamydomonas smithii)
UP000002008 324602  CHLAA      3850       0    3853  Chloroflexus aurantiacus (strain ATCC 29366 / DSM 635 / J-10-fl)
UP000008144 7719    CIOIN     16678     631   17317  Ciona intestinalis (Transparent sea squirt) (Ascidia intestinalis)
UP000002149 214684  CRYNJ      6603     143    6783  Cryptococcus neoformans var. neoformans serotype D (strain JEC21 / ATCC MYA-565) (Filobasidiella neoformans)
UP000000437 7955    DANRE     25043   18715   44479  Danio rerio (Zebrafish) (Brachydanio rerio)
UP000002524 243230  DEIRA      3085       0    3103  Deinococcus radiodurans (strain ATCC 13939 / DSM 20539 / JCM 16871 / LMG 4051 / NBRC 15346 / NCIMB 9279 / R1 / VKM B-1422)
UP000007719 515635  DICTD      1743       0    1744  Dictyoglomus turgidum (strain Z-1310 / DSM 6724)
UP000002195 44689   DICDI     12735      28   13387  Dictyostelium discoideum (Slime mold)
UP000000803 7227    DROME     13640    9651   23488  Drosophila melanogaster (Fruit fly)
UP000000625 83333   ECOLI      4306       9    8436  Escherichia coli (strain K12)
UP000002521 190304  FUSNN      2046       0    2064  Fusobacterium nucleatum subsp. nucleatum (strain ATCC 25586 / CIP 101130 / JCM 8532 / LMG 13131)
UP000000539 9031    CHICK     18557   11177   30403  Gallus gallus (Chicken)
UP000000577 243231  GEOSL      3402       0    3433  Geobacter sulfurreducens (strain ATCC 51573 / DSM 12127 / PCA)
UP000001548 184922  GIAIC      7154       0    7364  Giardia intestinalis (strain ATCC 50803 / WB clone C6) (Giardia lamblia)
UP000000557 251221  GLOVI      4406       0    4427  Gloeobacter violaceus (strain PCC 7421)
UP000001519 9595    GORGO     20946    6351   27477  Gorilla gorilla gorilla (Western lowland gorilla)
UP000000554 64091   HALSA      2426       0    2607  Halobacterium salinarum (strain ATCC 700922 / JCM 11081 / NRC-1) (Halobacterium halobium)
UP000000429 85962   HELPY      1553       0    1577  Helicobacter pylori (strain ATCC 700392 / 26695) (Campylobacter pylori)
UP000015101 6412    HELRO     23328       0   23433  Helobdella robusta (Californian leech)
UP000005640 9606    HUMAN     21042   71889   93380  Homo sapiens (Human)
UP000001555 6945    IXOSC     20469       4   20884  Ixodes scapularis (Black-legged tick) (Deer tick)
UP000001686 374847  KORCO      1602       0    1602  Korarchaeum cryptofilum (strain OPF8)
UP000000542 5664    LEIMA      8038       0    8354  Leishmania major
UP000018468 7918    LEPOC     18314    4146   22480  Lepisosteus oculatus (Spotted gar)
UP000001408 189518  LEPIN      3676       0    3707  Leptospira interrogans serogroup Icterohaemorrhagiae serovar Lai (strain 56601)
UP000000805 243232  METJA      1787       0    1788  Methanocaldococcus jannaschii (strain ATCC 43067 / DSM 2661 / JAL-1 / JCM 10045 / NBRC 100440) (Methanococcus jannaschii)
UP000002487 188937  METAC      4468       0    4540  Methanosarcina acetivorans (strain ATCC 35395 / DSM 2834 / JCM 12185 / C2A)
UP000002280 13616   MONDO     21271     973   22320  Monodelphis domestica (Gray short-tailed opossum)
UP000001357 81824   MONBE      9188       0    9204  Monosiga brevicollis (Choanoflagellate)
UP000000589 10090   MOUSE     22262   36810   59533  Mus musculus (Mouse)
UP000001584 83332   MYCTU      3993       4    5634  Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)
UP000000807 243273  MYCGE       483       0     483  Mycoplasma genitalium (strain ATCC 33530 / G-37 / NCTC 10195)
UP000000425 122586  NEIMB      2001       0    2062  Neisseria meningitidis serogroup B (strain MC58)
UP000001593 45351   NEMVE     24428      16   24785  Nematostella vectensis (Starlet sea anemone)
UP000002530 330879  ASPFU      9648       0    9653  Neosartorya fumigata (strain ATCC MYA-4609 / Af293 / CBS 101355 / FGSC A1100) (Aspergillus fumigatus)
UP000001805 367110  NEUCR      9759     508   17903  Neurospora crassa (strain ATCC 24698 / 74-OR23-1A / CBS 708.71 / DSM 1257 / FGSC 987)
UP000000792 436308  NITMS      1795       0    1795  Nitrosopumilus maritimus (strain SCM1)
UP000059680 39947   ORYSJ     44321    4922   52685  Oryza sativa subsp. japonica (Rice)
UP000001038 8090    ORYLA     19663    4982   24709  Oryzias latipes (Japanese rice fish) (Japanese killifish)
UP000002277 9598    PANTR     18980    1181   20273  Pan troglodytes (Chimpanzee)
UP000000600 5888    PARTE     39461       0   39753  Paramecium tetraurelia
UP000001055 321614  PHANO     15998       0   16008  Phaeosphaeria nodorum (strain SN15 / ATCC MYA-4574 / FGSC 10173) (Glume blotch fungus) (Parastagonospora nodorum)
UP000006727 3218    PHYPA     34813      24   35187  Physcomitrella patens subsp. patens (Moss)
UP000005238 164328  PHYRM     15349       0   15605  Phytophthora ramorum (Sudden oak death agent)
UP000001450 36329   PLAF7      5360       9    8104  Plasmodium falciparum (isolate 3D7)
UP000002438 208964  PSEAE      5562       2    5576  Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)
UP000008783 418459  PUCGT     15688     120   15920  Puccinia graminis f. sp. tritici (strain CRL 75-36-700-3 / race SCCL) (Black stem rust fungus)
UP000002494 10116   RAT       21412    9977   32826  Rattus norvegicus (Rat)
UP000001025 243090  RHOBA      7271       0    7325  Rhodopirellula baltica (strain DSM 10527 / NCIMB 13988 / SH1)
UP000002311 559292  YEAST      6722     122    6852  Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)
UP000002485 284812  SCHPO      5142       9    5155  Schizosaccharomyces pombe (strain 972 / ATCC 24843) (Fission yeast)
UP000001312 665079  SCLS1     14445       0   16092  Sclerotinia sclerotiorum (strain ATCC 18683 / 1980 / Ss-1) (White mold) (Whetzelinia sclerotiorum)
UP000001973 100226  STRCO      8038       1    8154  Streptomyces coelicolor (strain ATCC BAA-471 / A3(2) / M145)
UP000001974 273057  SULSO      2938       0    2995  Sulfolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2)
UP000001425 1111708 SYNY3      3507       0    3566  Synechocystis sp. (strain PCC 6803 / Kazusa)
UP000001449 35128   THAPS     11717       1   23292  Thalassiosira pseudonana (Marine diatom) (Cyclotella nana)
UP000000536 69014   THEKO      2301       0    2305  Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1))
UP000000718 289376  THEYD      1982       0    2033  Thermodesulfovibrio yellowstonii (strain ATCC 51303 / DSM 11347 / YP87)
UP000008183 243274  THEMA      1852       0    2872  Thermotoga maritima (strain ATCC 43589 / MSB8 / DSM 3109 / JCM 10099)
UP000007266 7070    TRICA     16563    1942   18563  Tribolium castaneum (Red flour beetle)
UP000001542 5722    TRIVA     50190       0   59680  Trichomonas vaginalis
UP000000561 237631  USTMA      6788      18    6833  Ustilago maydis (strain 521 / FGSC 9021) (Corn smut fungus)
UP000008143 8364    XENTR     24177    5510   29854  Xenopus tropicalis (Western clawed frog) (Silurana tropicalis)
UP000001300 284591  YARLI      6448       5    6476  Yarrowia lipolytica (strain CLIB 122 / E 150) (Yeast) (Candida lipolytica)
UP000007305 4577    MAIZE     39476   60000  100392  Zea mays (Maize)


Gene mapping files (*.gene2acc)
===============================

Column 1 is a unique gene symbol that is chosen with the following order of
preference from the annotation found in:
1) Model Organism Database (MOD)
2) Ensembl or Ensembl Genomes database
3) UniProt Ordered Locus Name (OLN)
4) UniProt Open Reading Frame (ORF)
5) UniProt Gene Name
A dash symbol ('-') is used when the gene encoding a protein is unknown.

Column 2 is the UniProtKB accession or isoform identifier for the given gene
symbol. This column may have redundancy when two or more genes have identical
translations.


Protein FASTA files (*.fasta and *_additional.fasta)
====================================================

These files, composed of canonical and additional sequences, are non-redundant
FASTA sets for the sequences of each reference proteome.
The additional set contains isoform/variant sequences for a given gene, and its
FASTA header indicates the corresponding canonical sequence ("Isoform of ...").
The FASTA format is the standard UniProtKB format.

For further references about the standard UniProtKB format, please see:
    http://www.uniprot.org/help/fasta-headers
    http://www.uniprot.org/faq/38

E.g. Canonical set:
>sp|Q9H6Y5|MAGIX_HUMAN PDZ domain-containing protein MAGIX OS=Homo sapiens GN=MAGIX PE=1 SV=3
MEPRTGGAANPKGSRGSRGPSPLAGPSARQLLARLDARPLAARAAVDVAALVRRAGATLR
LRRKEAVSVLDSADIEVTDSRLPHATIVDHRPQHRWLETCNAPPQLIQGKAHSAPKPSQA
SGHFSVELVRGYAGFGLTLGGGRDVAGDTPLAVRGLLKDGPAQRCGRLEVGDVVLHINGE
STQGLTHAQAVERIRAGGPQLHLVIRRPLETHPGKPRGVGEPRKGVVPSWPDRSPDPGGP
EVTGSRSSSTSLVQHPPSRTTLKKTRGSPEPSPEAAADGPTVSPPERRAEDPNDQIPGSP
GPWLVPSEERLSRALGVRGAAQFAQEMAAGRRRH

E.g. Additional sets:
>sp|Q9H6Y5-2|MAGIX_HUMAN Isoform of Q9H6Y5, Isoform 2 of PDZ domain-containing protein MAGIX OS=Homo sapiens GN=MAGIX
MPLLWITGPRYHLILLSEASCLRANYVHLCPLFQHRWLETCNAPPQLIQGKAHSAPKPSQ
ASGHFSVELVRGYAGFGLTLGGGRDVAGDTPLAVRGLLKDGPAQRCGRLEVGDVVLHING
ESTQGLTHAQAVERIRAGGPQLHLVIRRPLETHPGKPRGVGEPRKGVVPSWPDRSPDPGG
PEVTGSRSSSTSLVQHPPSRTTLKKTRGSPEPSPEAAADGPTVSPPERRAEDPNDQIPGS
PGPWLVPSEERLSRALGVRGAAQFAQEMAAGRRRH
>tr|C9J123|C9J123_HUMAN Isoform of Q9H6Y5, PDZ domain-containing protein MAGIX (Fragment) OS=Homo sapiens GN=MAGIX PE=1 SV=2
MSPNSPLHCFYLPAVSVLDSADIEVTDSRLPHATIVDHRPQVGDLVLHINGESTQGLTHA
QAVERIRAGGPQLHLVIRRPLETHPGKPRGVGEPRKGVDRSPDPGGPEVTGSRSSSTSLV
QHPPSRTTLKKTRGSPEPSPEAA


Coding DNA Sequence FASTA files (*_DNA.fasta)
=============================================

These files contain the coding DNA sequences (CDS) for the protein sequences
where it was possible.
The format is as in the following example (UP000005640_9606_DNA.fasta):

>sp|A0A183|ENSP00000411070
ATGTCACAGCAGAAGCAGCAATCTTGGAAGCCTCCAAATGTTCCCAAATGCTCCCCTCCC
CAAAGATCAAACCCCTGCCTAGCTCCCTACTCGACTCCTTGTGGTGCTCCCCATTCAGAA
GGTTGTCATTCCAGTTCCCAAAGGCCTGAGGTTCAGAAGCCTAGGAGGGCTCGTCAAAAG
CTGCGCTGCCTAAGTAGGGGCACAACCTACCACTGCAAAGAGGAAGAGTGTGAAGGCGAC
TGA

The 3 fields of the FASTA header are:
1) 'sp' (Swiss-Prot reviewed) or 'tr' (TrEMBL)
2) UniProtKB Accession
3) EMBL Protein ID or Ensembl/Ensembl Genome ID


Unsuccessful Coding DNA Sequence mapping files (*_DNA.miss)
===========================================================

For the species that did not have a perfect mapping for all protein sequences
to a CDS, these files contain the entries that could not be mapped.
The format is as in the following example (UP000005640_9606_DNA.miss):

sp      A6NF01  CAUTION: Could be the product of a pseudogene.
sp      A6NFI3  NOT_ANNOTATED_CDS

The 3 fields are:
1) 'sp' (Swiss-Prot reviewed) or 'tr' (TrEMBL)
2) UniProtKB accession
3) Reason why the protein could not be mapped to a CDS


Database mapping files (*.idmapping)
====================================

These files contain mappings from UniProtKB to other databases for each
reference proteome.
The format consists of three tab-separated columns:

1. UniProtKB accession
2. ID_type:
   Database name as shown in UniProtKB cross-references and supported by the ID
   mapping tool on the UniProt web site (http://www.uniprot.org/mapping)
3. ID:
   Identifier in the cross-referenced database.


SeqXML files (*.xml)
====================

The xml files contain all the information from fasta (canonical and additional),
idmapping and CDS in SeqXML format (see http://seqxml.org.)
E.g. (from UP000005640_9606.xml, header and one entry)

<?xml version="1.0" encoding="utf-8"?>
<seqXML xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" speciesName="Homo sapiens" xsi:noNamespaceSchemaLocation="http://www.seqxml.org/0.4/seqxml.xsd" seqXMLversion="0.4" sourceVersion="2016_04" source="QfO http://w
ww.ebi.ac.uk/reference_proteomes/" ncbiTaxID="9606">
  <entry source="UniProtKB" id="U3KQE9">
    <description>Uncharacterized protein (Fragment)</description>
    <AAseq>XLLLAINGVTECFTFAAMSKEEVDRYNFVMLALSSSFLVLSYLLTRWCGSVGFILANCFNMGIRITQSLCFIHRYYRRSPHRPLAGLHLSPVLLGTFALSGGVTAVSERPRKPWRSSGPWCPC</AAseq>
    <DBRef source="UniProtKB-ID" id="U3KQE9_HUMAN"/>
    <DBRef source="UniRef100" id="UniRef100_U3KQE9"/>
    <DBRef source="UniRef90" id="UniRef90_U3KQE9"/>
    <DBRef source="UniRef50" id="UniRef50_U3KQE9"/>
    <DBRef source="UniParc" id="UPI00038BAF87"/>
    <DBRef source="EMBL" id="AC096887"/>
    <DBRef source="EMBL" id="AC099667"/>
    <DBRef source="EMBL" id="AC103589"/>
    <DBRef source="EMBL-CDS" id="-"/>
    <DBRef source="NCBI_TaxID" id="9606"/>
    <DBRef source="STRING" id="9606.ENSP00000296292"/>
    <DBRef source="Ensembl" id="ENSG00000272305"/>
    <DBRef source="Ensembl_TRS" id="ENST00000607283"/>
    <DBRef source="Ensembl_PRO" id="ENSP00000475819"/>
    <DBRef source="UCSC" id="uc062krn.1"/>
    <DBRef source="eggNOG" id="KOG2864"/>
    <DBRef source="eggNOG" id="ENOG410Y1D5"/>
    <DBRef source="GeneTree" id="ENSGT00390000011390"/>
    <DBRef source="OMA" id="ERPRKPW"/>
    <DBRef source="UniProt" id="U3KQG8"/>
    <DBRef source="UniProt" id="U3KQK6"/>
    <property name="DNAsource" value="ENSP00000475819"/>
    <property name="ensemblVersion" value="84,31"/>
    <property name="UPID" value="UP000005640"/>
    <property name="SV" value="1"/>
    <property name="DNAseq" value="NTTCTCCTGCTTGCCATCAATGGAGTGACAGAGTGTTTCACATTTGCTGCCATGAGCAAAGAGGAGGTCGACAGGTACAATTTTGTGATGCTGGCCCTGTCCTCCTCATTCCTGGTGTTATCCTATCTCTTGACCCGTTGGTGTGGCAGCGTGGGCTTCATCTTGGCCAACTGCTTTAACATGGGCATTCGGATCACGCAGAGCCTTTGCTTCATCCACCGCTACTACCGAAGGAGCCCCCACAGGCCCCTGGCTGGCCTGCACCTATCGCCAGTCCTGCTCGGGACATTTGCCCTCAGTGGTGGGGTTACTGCTGTTTCGGAGAGGCCAAGAAAACCATGGAGGAGCAGTGGACCTTGGTGTCCCTGCTGA"/>
    <property name="PE" value="4"/>
  </entry>


--------------------------------------------------------------------------------
  LICENSE
--------------------------------------------------------------------------------
We have chosen to apply the Creative Commons Attribution-NoDerivs License to all
copyrightable parts of our databases. This means that you are free to copy,
distribute, display and make commercial use of these databases in all
legislations, provided you give us credit. However, if you intend to distribute
a modified version of one of our databases, you must ask us for permission
first.

(c) 2002-2017 UniProt Consortium

--------------------------------------------------------------------------------
  DISCLAIMER
--------------------------------------------------------------------------------
We make no warranties regarding the correctness of the data, and disclaim
liability for damages resulting from its use. We cannot provide unrestricted
permission regarding the use of the data, as some data may be covered by patents
or other rights.

Any medical or genetic information is provided for research, educational and
informational purposes only. It is not in any way intended to be used as a
substitute for professional medical advice, diagnosis, treatment or care.

### Included Genomes ###
UP000007062 7165 ANOGA 11930 1442 14384 Anopheles gambiae (African malaria mosquito)
UP000000798 224324 AQUAE 1553 0 1557 Aquifex aeolicus (strain VF5)
UP000006548 3702 ARATH 27502 13831 41673 Arabidopsis thaliana (Mouse-ear cress)
UP000001570 224308 BACSU 4197 7 4205 Bacillus subtilis (strain 168)
UP000001414 226186 BACTN 4782 0 4823 Bacteroides thetaiotaomicron (strain ATCC 29148 / DSM 2079 / NCTC 10582 / E50 / VPI-5482)
UP000007241 684364 BATDJ 8610 0 8700 Batrachochytrium dendrobatidis (strain JAM81 / FGSC 10211) (Frog chytrid fungus)
UP000009136 9913 BOVIN 21987 2529 24612 Bos taurus (Bovine)
UP000002526 224911 BRADU 8253 0 8314 Bradyrhizobium diazoefficiens (strain JCM 10833 / IAM 13628 / NBRC 14792 / USDA 110)
UP000001554 7739 BRAFL 28542 2 28639 Branchiostoma floridae (Florida lancelet) (Amphioxus)
UP000001940 6239 CAEEL 20057 7998 28264 Caenorhabditis elegans
UP000000559 237561 CANAL 6153 1 10852 Candida albicans (strain SC5314 / ATCC MYA-2876) (Yeast)
UP000002254 9615 CANLF 20141 5395 25631 Canis lupus familiaris (Dog) (Canis familiaris)
UP000000431 272561 CHLTR 895 0 895 Chlamydia trachomatis (strain D/UW-3/Cx)
UP000006906 3055 CHLRE 14271 67 14483 Chlamydomonas reinhardtii (Chlamydomonas smithii)
UP000002008 324602 CHLAA 3850 0 3853 Chloroflexus aurantiacus (strain ATCC 29366 / DSM 635 / J-10-fl)
UP000008144 7719 CIOIN 16678 631 17317 Ciona intestinalis (Transparent sea squirt) (Ascidia intestinalis)
UP000002149 214684 CRYNJ 6603 143 6783 Cryptococcus neoformans var. neoformans serotype D (strain JEC21 / ATCC MYA-565) (Filobasidiella neoformans)
UP000000437 7955 DANRE 25043 18715 44479 Danio rerio (Zebrafish) (Brachydanio rerio)
UP000002524 243230 DEIRA 3085 0 3103 Deinococcus radiodurans (strain ATCC 13939 / DSM 20539 / JCM 16871 / LMG 4051 / NBRC 15346 / NCIMB 9279 / R1 / VKM B-1422)
UP000007719 515635 DICTD 1743 0 1744 Dictyoglomus turgidum (strain Z-1310 / DSM 6724)
UP000002195 44689 DICDI 12735 28 13387 Dictyostelium discoideum (Slime mold)
UP000000803 7227 DROME 13640 9651 23488 Drosophila melanogaster (Fruit fly)
UP000000625 83333 ECOLI 4306 9 8436 Escherichia coli (strain K12)
UP000002521 190304 FUSNN 2046 0 2064 Fusobacterium nucleatum subsp. nucleatum (strain ATCC 25586 / CIP 101130 / JCM 8532 / LMG 13131)
UP000000539 9031 CHICK 18557 11177 30403 Gallus gallus (Chicken)
UP000000577 243231 GEOSL 3402 0 3433 Geobacter sulfurreducens (strain ATCC 51573 / DSM 12127 / PCA)
UP000001548 184922 GIAIC 7154 0 7364 Giardia intestinalis (strain ATCC 50803 / WB clone C6) (Giardia lamblia)
UP000000557 251221 GLOVI 4406 0 4427 Gloeobacter violaceus (strain PCC 7421)
UP000001519 9595 GORGO 20946 6351 27477 Gorilla gorilla gorilla (Western lowland gorilla)
UP000000554 64091 HALSA 2426 0 2607 Halobacterium salinarum (strain ATCC 700922 / JCM 11081 / NRC-1) (Halobacterium halobium)
UP000000429 85962 HELPY 1553 0 1577 Helicobacter pylori (strain ATCC 700392 / 26695) (Campylobacter pylori)
UP000015101 6412 HELRO 23328 0 23433 Helobdella robusta (Californian leech)
UP000005640 9606 HUMAN 21042 71889 93380 Homo sapiens (Human)
UP000001555 6945 IXOSC 20469 4 20884 Ixodes scapularis (Black-legged tick) (Deer tick)
UP000001686 374847 KORCO 1602 0 1602 Korarchaeum cryptofilum (strain OPF8)
UP000000542 5664 LEIMA 8038 0 8354 Leishmania major
UP000018468 7918 LEPOC 18314 4146 22480 Lepisosteus oculatus (Spotted gar)
UP000001408 189518 LEPIN 3676 0 3707 Leptospira interrogans serogroup Icterohaemorrhagiae serovar Lai (strain 56601)
UP000000805 243232 METJA 1787 0 1788 Methanocaldococcus jannaschii (strain ATCC 43067 / DSM 2661 / JAL-1 / JCM 10045 / NBRC 100440) (Methanococcus jannaschii)
UP000002487 188937 METAC 4468 0 4540 Methanosarcina acetivorans (strain ATCC 35395 / DSM 2834 / JCM 12185 / C2A)
UP000002280 13616 MONDO 21271 973 22320 Monodelphis domestica (Gray short-tailed opossum)
UP000001357 81824 MONBE 9188 0 9204 Monosiga brevicollis (Choanoflagellate)
UP000000589 10090 MOUSE 22262 36810 59533 Mus musculus (Mouse)
UP000001584 83332 MYCTU 3993 4 5634 Mycobacterium tuberculosis (strain ATCC 25618 / H37Rv)
UP000000807 243273 MYCGE 483 0 483 Mycoplasma genitalium (strain ATCC 33530 / G-37 / NCTC 10195)
UP000000425 122586 NEIMB 2001 0 2062 Neisseria meningitidis serogroup B (strain MC58)
UP000001593 45351 NEMVE 24428 16 24785 Nematostella vectensis (Starlet sea anemone)
UP000002530 330879 ASPFU 9648 0 9653 Neosartorya fumigata (strain ATCC MYA-4609 / Af293 / CBS 101355 / FGSC A1100) (Aspergillus fumigatus)
UP000001805 367110 NEUCR 9759 508 17903 Neurospora crassa (strain ATCC 24698 / 74-OR23-1A / CBS 708.71 / DSM 1257 / FGSC 987)
UP000000792 436308 NITMS 1795 0 1795 Nitrosopumilus maritimus (strain SCM1)
UP000059680 39947 ORYSJ 44321 4922 52685 Oryza sativa subsp. japonica (Rice)
UP000001038 8090 ORYLA 19663 4982 24709 Oryzias latipes (Japanese rice fish) (Japanese killifish)
UP000002277 9598 PANTR 18980 1181 20273 Pan troglodytes (Chimpanzee)
UP000000600 5888 PARTE 39461 0 39753 Paramecium tetraurelia
UP000001055 321614 PHANO 15998 0 16008 Phaeosphaeria nodorum (strain SN15 / ATCC MYA-4574 / FGSC 10173) (Glume blotch fungus) (Parastagonospora nodorum)
UP000006727 3218 PHYPA 34813 24 35187 Physcomitrella patens subsp. patens (Moss)
UP000005238 164328 PHYRM 15349 0 15605 Phytophthora ramorum (Sudden oak death agent)
UP000001450 36329 PLAF7 5360 9 8104 Plasmodium falciparum (isolate 3D7)
UP000002438 208964 PSEAE 5562 2 5576 Pseudomonas aeruginosa (strain ATCC 15692 / DSM 22644 / CIP 104116 / JCM 14847 / LMG 12228 / 1C / PRS 101 / PAO1)
UP000008783 418459 PUCGT 15688 120 15920 Puccinia graminis f. sp. tritici (strain CRL 75-36-700-3 / race SCCL) (Black stem rust fungus)
UP000002494 10116 RAT 21412 9977 32826 Rattus norvegicus (Rat)
UP000001025 243090 RHOBA 7271 0 7325 Rhodopirellula baltica (strain DSM 10527 / NCIMB 13988 / SH1)
UP000002311 559292 YEAST 6722 122 6852 Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)
UP000002485 284812 SCHPO 5142 9 5155 Schizosaccharomyces pombe (strain 972 / ATCC 24843) (Fission yeast)
UP000001312 665079 SCLS1 14445 0 16092 Sclerotinia sclerotiorum (strain ATCC 18683 / 1980 / Ss-1) (White mold) (Whetzelinia sclerotiorum)
UP000001973 100226 STRCO 8038 1 8154 Streptomyces coelicolor (strain ATCC BAA-471 / A3(2) / M145)
UP000001974 273057 SULSO 2938 0 2995 Sulfolobus solfataricus (strain ATCC 35092 / DSM 1617 / JCM 11322 / P2)
UP000001425 1111708 SYNY3 3507 0 3566 Synechocystis sp. (strain PCC 6803 / Kazusa)
UP000001449 35128 THAPS 11717 1 23292 Thalassiosira pseudonana (Marine diatom) (Cyclotella nana)
UP000000536 69014 THEKO 2301 0 2305 Thermococcus kodakarensis (strain ATCC BAA-918 / JCM 12380 / KOD1) (Pyrococcus kodakaraensis (strain KOD1))
UP000000718 289376 THEYD 1982 0 2033 Thermodesulfovibrio yellowstonii (strain ATCC 51303 / DSM 11347 / YP87)
UP000008183 243274 THEMA 1852 0 2872 Thermotoga maritima (strain ATCC 43589 / MSB8 / DSM 3109 / JCM 10099)
UP000007266 7070 TRICA 16563 1942 18563 Tribolium castaneum (Red flour beetle)
UP000001542 5722 TRIVA 50190 0 59680 Trichomonas vaginalis
UP000000561 237631 USTMA 6788 18 6833 Ustilago maydis (strain 521 / FGSC 9021) (Corn smut fungus)
UP000008143 8364 XENTR 24177 5510 29854 Xenopus tropicalis (Western clawed frog) (Silurana tropicalis)
UP000001300 284591 YARLI 6448 5 6476 Yarrowia lipolytica (strain CLIB 122 / E 150) (Yeast) (Candida lipolytica)
UP000007305 4577 MAIZE 39476 60000 100392 Zea mays (Maize)
