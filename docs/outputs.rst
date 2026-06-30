Outputs
========

Annotation
-----------

Arraystar
^^^^^^^^^

.. code-block:: text

    arraystar       hg19                            hg38                            circRNA Alias           circbase              source          chrom   strand  txStart         txEnd         circRNA_type    best_transcript GeneSymbol
    ASCRP3010903    chr1:101372407-101387397|+      chr1:100906851-100921841|+      hsa_circRNA_100291      hsa_circ_0000098      circBase        chr1    +       101372407       101387397     exonic          NM_133496       SLC30A7
    ASCRP3002028    chr1:11130956-11137005|-        chr1:11070899-11076948|-        hsa_circRNA_100050      hsa_circ_0009759      circBase        chr1    -       11130956        11137005      exonic          NM_002685       EXOSC10
    ASCRP3004170    chr1:113143415-113153625|-      chr1:112600793-112611003|-      hsa_circRNA_000108      hsa_circ_0000108      circBase        chr1    -       113143415       113153625     exonic          NM_017744       ST7L

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - arraystar
     - The Arraystar circRNA ID.
   * - hg19
     - The circRNA coordinates in the hg19 reference genome.
   * - hg38
     - The circRNA coordinates in the hg38 reference genome.
   * - Alias
     - The circRNA alias. (does not align with other databases)
   * - circbase
     - circBase identifier for the circRNA.
   * - source
     - Database source of the circRNA annotation.
   * - chrom
     - Chromosome of the circRNA.
   * - strand
     - Strand of the circRNA.
   * - txStart
     - Start coordinate of the circRNA.
   * - txEnd
     - End coordinate of the circRNA.
   * - circRNA_type
     - Type of the circRNA (antisense, exonic, intergenic, intronic, sense overlapping).
   * - best_transcript
     - transcript ID from which the circRNA likely originates. Null if the circRNA is intergenic.
   * - GeneSymbol
     - Gene symbol of the overlapping gene for the circRNA. Null if the circRNA is intergenic.


CircAtlas
^^^^^^^^^

Please note that circatlas annotations are not readily available publicly. For this reason, mapping keys between circatlas, hg19 and hg38 are provided as outputs:

.. code-block:: text

    circatlas       hg19                    hg38
    hsa-SDF4_0001   chr1:1158623-1159348|-  chr1:1223243-1223968|-
    hsa-SKI_0001    chr1:2234416-2236024|+  chr1:2302977-2304585|+
    hsa-RERE_0008   chr1:8555122-8674745|-  chr1:8495062-8614686|-


Circbank
^^^^^^^^

.. code-block:: text

    circbank                hg19                            hg38                            aliasID                                                     exonNumber     Exons                                                           circType        splicedLength   Assembly      geneID                  geneSymbol      bestTranscriptID        Biotype 
    hsa_SLC30A7_0001800     chr1:101372407-101387397|+      chr1:100906851-100921841|+      hsa_circ_0000098,hsa_circSLC30A7_010,hsa-SLC30A7_0001       6.0            101372407-101372521,101376618-101376706,101379218-101379362     Exon-Exon       660             hg19          ENSG00000162695.7       SLC30A7         ENST00000370112.4       protein_coding
    hsa_EXOSC10_0002000     chr1:11130956-11137005|-        chr1:11070899-11076948|-        hsa_circ_0009759,hsa_circEXOSC10_023,hsa-EXOSC10_0001       5.0            11130956-11131030,11132143-11132228,11133990-11134065,          Exon-Exon       437             hg19          ENSG00000171824.9       EXOSC10         ENST00000376936.4       protein_coding
    hsa_ST7L_0015700        chr1:113143415-113153625|-      chr1:112600793-112611003|-      hsa_circ_0000108,hsa_circST7L_014,hsa-ST7L_0039             2.0            113143415-113143470,113153462-113153625                         Exon-Exon       218             hg19          ENSG00000007341.14      ST7L            ENST00000369664.1       protein_coding



.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - circbank
     - The circBank circRNA ID.
   * - hg19
     - The circRNA coordinates in the hg19 reference genome.
   * - hg38
     - The circRNA coordinates in the hg38 reference genome.
   * - aliasID
     - Three aliases: circBank ID, CircBank ID (alias), and circAtlas ID.
   * - exonNumber
     - Number of exons in the circRNA.
   * - Exons
     - Comma separated list of exon start-end coordinates for the circRNA.
   * - circType
     - Type of the circRNA (Intergenic, Intron, Exonic or any combination thereof)
   * - splicedLength
     - Length of the spliced circRNA.
   * - Assembly
      - Reference genome assembly used by circbank.
   * - geneID
      - Ensembl gene version ID of the overlapping gene for the circRNA.
   * - geneSymbol
      - Gene symbol of the overlapping gene for the circRNA.
   * - bestTranscriptID
      - Ensembl transcript version ID from which the circRNA likely originates.
   * - Biotype
      - Biotype of the best transcript. Null if the circRNA is intergenic.

 


circBase
^^^^^^^^

.. code-block:: text
    
    circbase              hg19                      hg38                      chrom    chromStart      chromEnd        score   strand  thickStart      thickEnd     itemRGB    blockCount      blockSizes         blockStarts
    hsa_circ_0000002      chr1:1158623-1159348|-    chr1:1223243-1223968|-    chr1     1158623         1159348         1000    -       1158623         1159348      0,0,255    2               114,137            0,588
    hsa_circ_0009244      chr1:1158623-1164326|-    chr1:1223243-1228946|-    chr1     1158623         1164326         1000    -       1158623         1164173      0,0,255    3               114,137,479        0,588,5224
    hsa_circ_0007120      chr1:2234416-2236024|+    chr1:2302977-2304585|+    chr1     2234416         2236024         1000    +       2234416         2236024      255,0,0    4               126,116,263,293    0,307,862,1315



.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
    * - circbase
     - The circBase circRNA ID.
    * - hg19
     - The circRNA coordinates in the hg19 reference genome.
    * - hg38
     - The circRNA coordinates in the hg38 reference genome.
    * - chrom
     - Chromosome of the circRNA.
    * - chromStart
     - Start coordinate of the circRNA.
    * - chromEnd
     - End coordinate of the circRNA.
    * - score
     - Score of the circRNA (not used).
    * - strand
     - Strand of the circRNA.
    * - thickStart
     - Start coordinate of the circRNA (same as chromStart).
    * - thickEnd
     - End coordinate of the circRNA (same as chromEnd).
    * - itemRGB
     - RGB color code for the circRNA (not used).
    * - blockCount
     - Number of exons in the circRNA.
    * - blockSizes
     - Comma separated list of exon sizes for the circRNA.
    * - blockStarts
     - Comma separated list of exon start coordinates relative to chromStart for the circRNA.


circBase outputs can be reverted to their original BED12 format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1):

.. code-block:: bash

    awk 'BEGIN{FS=OFS="\t"} {
      out=$4 OFS $5 OFS $6 OFS $1
      for(i=7;i<=NF;i++) out=out OFS $i
      print out
    }' circbase_hits.txt > circbase_bed12.bed


CIRCpedia
^^^^^^^^^

.. code-block:: text

    circpedia           hg19                          hg38                          circname        gene_Ensembl            gene_Refseq     transcript_Ensembl      transcript_Refseq
    CIRCHSA_PUM1_99     chr1:31465236-31468067|-      chr1:30992389-30995220|-      circPUM1(6,7)   ENSG00000134644.16      PUM1            ENST00000257075.9       PUM1-201
    CIRCHSA_SUCLG1_4    chr2:84652538-84676876|-      chr2:84425414-84449752|-      circSUCLG1(2-8) ENSG00000163541.12      SUCLG1          ENST00000393868.7       SUCLG1-201
    CIRCHSA_SOGA1_5     chr20:35421651-35467844|-     chr20:36793248-36839441|-     circSOGA1(2-14) ENSG00000149639.15      SOGA1           ENST00000237536.9       SOGA1-201



.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - circpedia
     - The CIRCpedia circRNA ID.
   * - hg19
     - The circRNA coordinates in the hg19 reference genome.
   * - hg38
     - The circRNA coordinates in the hg38 reference genome.
   * - circname
     - CIRCpedia naming convention (see `PMC10114414 <https://pmc.ncbi.nlm.nih.gov/articles/PMC10114414/>`_)
   * - gene_Ensembl
     - Ensembl gene version ID of the overlapping gene for the circRNA.
   * - gene_Refseq
     - RefSeq gene ID of the overlapping gene for the circRNA.
   * - transcript_Ensembl
     - Ensembl transcript version ID from which the circRNA likely originates.
   * - transcript_Refseq
     - RefSeq transcript ID from which the circRNA likely originates.


circRNADb
^^^^^^^^^

.. code-block:: text

  circRNA_DB      hg19    hg38    gene_symbol     genomic_length  transcript_id   spliced_length  exon_number     exon_sizes      exon_offsets    samples pubmed_id
  hsa_circ_00089  chr13:42385360-42393522|-       chr13:41811224-41819386|-       VWA8    8162    NM_001009814    363     3       116,78,169      0,5473,7993     H9 hESCs,normal brain tissue,oligodendroma,Hs68  25242744,26873924,23249747
  hsa_circ_00131  chr16:69189773-69201088|+       chr16:69155870-69167185|+       CIRH1A  11315   NM_032830       780     6       123,157,107,96,186,111  0,1213,4485,7212,9470,11204     H9 hESCs,Hs68,leukemia   25242744,23249747,22319583
  hsa_circ_00164  chr2:20507738-20527139|-        chr2:20307977-20327378|-        PUM2    19401   NM_015317       901     6       94,271,170,188,109,69   0,336,3516,4258,10559,19332     H9 hESCs 25242744

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - circRNA_DB
     - The circRNADb circRNA ID, not to be confused with circBase IDs.
   * - hg19
     - The circRNA coordinates in the hg19 reference genome.
   * - hg38
     - The circRNA coordinates in the hg38 reference genome.
   * - gene_symbol
     - Gene symbol of the overlapping gene for the circRNA.
   * - genomic_length
     - Length of the circRNA in the genome.
   * - transcript_id
     - RefSeq transcript ID from which the circRNA likely originates.
   * - spliced_length
     - Length of the spliced circRNA.
   * - exon_number
     - Number of exons in the circRNA.
   * - exon_sizes
     - Comma separated list of exon sizes for the circRNA.
   * - exon_offsets
     - Comma separated list of exon start coordinates relative to the circRNA start coordinate.
   * - samples
     - Comma separated list of samples in which the circRNA was detected.
   * - pubmed_id
     - Comma separated list of PubMed IDs for publications in which the circRNA was detected.

CSCD
^^^^

.. code-block:: text

  hg38                            hg19                                circbase_id             type            position        host_gene   read_counts   algorithm                                     subcellular_location                                                             gene_type   alternative_splicing    ratios          works_cited    mioncocirc      number_of_algorithm
  chr1:100906851-100921841|+      chr1:101372407-101387397|+          hsa_circ_0000098        normal,cancer   Exon,Intron     SLC30A7     SEE_BELOW     circRNA_finder,CIRI,find_circ,circexplorer2   nucleoplasmic,chromatin,cytosolic,nuclear,membrane,insoluble_cytoplasmic,---     mRNA        SEE_BELOW               3.2727267736    Marcel_2019    Mioncocirc      4
  chr1:108699570-108699721|+      chr1:109242192-109242343|+                                  normal,cancer   Exon            PRPF38B     SEE_BELOW     circRNA_finder,find_circ                      nucleoplasmic,nucleolus,chromatin,                                               mRNA        SEE_BELOW               1.1759283247                                   2
  chr1:100424221-100442996|+      chr1-100889777-100908552|+          hsa_circ_0000097        normal,cancer   Exon,Intron     CDC14A      SEE_BELOW     circRNA_finder,CIRI,find_circ,circexplorer2   nucleoplasmic,nucleolus,chromatin,cytosolic,nuclear,membrane,                    mRNA        SEE_BELOW               2.6381546935    Marcel_2019    Mioncocirc      4

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - hg38
     - The circRNA coordinates in the hg38 reference genome.
   * - hg19
     - The circRNA coordinates in the hg19 reference genome.
   * - circbase_id
     - circBase identifier for the circRNA if available.
   * - type
     - Type of the circRNA (normal, cancer, or both).
   * - position
     - Position of the circRNA (Exon, Intron, Intergenic or Exon,Intron i.e EIcircRNA).
   * - host_gene
     - Gene symbol of the overlapping gene for the circRNA.
   * - read_counts
     - Read counts for the circRNA; truncated in output, full explination below table.
   * - algorithm
     - Algorithms that detected the circRNA.
   * - subcellular_location
     - Subcellular locations of the circRNA.
   * - gene_type
     - Gene type of the overlapping gene for the circRNA.
   * - alternative_splicing
     - Alternative splcing events see below for full description.
   * - ratios
     - Ratios of circRNA to linear RNA expression.
   * - works_cited
     - Citation for the circRNA.
   * - mioncocirc
     - Presence of the circRNA in the Mioncocirc database.
   * - number_of_algorithm
     - Number of algorithms that detected the circRNA.

Alternative splicing events are reported in the following format:

.. code-block:: text

    <start coordinate>:<end coordinate>,<event type>;<start coordinate>:<end coordinate>,<event type>;...
    100436132:100436242,ES;100436132:100436242,ES;100431324:100431361,ES;

Where event type is one of the following:

.. list-table::
   :header-rows: 1
   :widths: 30 
  * - Event Type
    - Description
  * - ES
    - Exon skipping
  * - A5SS
    - Alternative 5' splice site
  * - A3SS
    - Alternative 3' splice site
  * - MXR
    - Mutually exclusive exons
  * - IR
    - Retained Intron 

Read counts are provided for circular RNAs in the following format:

.. code-block:: text

    <sample_name>,<algorithm>,<read_count> | <sample_name>,<algorithm>,<read_count> | ...
    mesenchymal_stem_cell_of_adipose_2,circexplorer2,1 | mesenchymal_stem_cell_of_adipose_2,circRNA_finder,1 | HCC_2,CIRI,2 | HCC_3,circexplorer2,2 | HCC_3,circRNA_finder,2 | HCC_3,CIRI,2

We recommend extracting the hg38/hg19 column and the read counts column into a long format dataframe:

.. list-table::
   :header-rows: 1
   :widths: 30 40 30 30 

   * - circRNA
     - Sample Name
     - Algorithm
     - Read Count
   * - chr1:100424221-100442996|+
     - mesenchymal_stem_cell_of_adipose_2
     - circexplorer2
     - 1
   * - chr1:100424221-100442996|+
     - mesenchymal_stem_cell_of_adipose_2
     - circRNA_finder
     - 1
   * - chr1:100906851-100921841|+
     - HCC_2
     - CIRI
     - 2
   * - chr1:100906851-100921841|+
     - HCC_3
     - circexplorer2
     - 2
   * - chr1:100906851-100921841|+
     - HCC_3
     - circRNA_finder
     - 2
   * - chr1:100906851-100921841|+
     - HCC_3
     - CIRI
     - 2


In this way, you can choose to generate a circRNA count matrix for common samples across algorithms by averaging their counts.

exoRbase
^^^^^^^^

FASTA
-----

Please note that the FASTA ouputs below are truncated for brevity.

Arraystar
^^^^^^^^^

.. code-block:: text

  >ASCRP3010903
  TAGTTGTAAGCTTAGGCTTGATTTCCGACTCTTTTCACATGTTTTTCGATAGCACTGCCA

Circbank
^^^^^^^^

.. code-block:: text

  >hsa_SLC30A7_0001800
  CTTAGGCTTGATTTCCGACTCTTTTCACATGTTTTTCGATAGCACTGCCATTTT

circBase
^^^^^^^^

.. code-block:: text

  >hsa_circ_0000002
  GGTGGATGTGAACACTGACCGGAAGATCAGTG

CIRCpedia
^^^^^^^^^

.. code-block:: text

  >CIRCHSA_SDF4_1
  GGTGGATGTGAACACTGACCGGAAGATCAG

circRNADb
^^^^^^^^^

.. code-block:: text

  >hsa_circ_00089
  ATCCATTTTTCCTATCCATC

CSCD
^^^^

.. code-block:: text

  >chr1:100424221-100442996|+
  CTTAGGCTTGATTTCCGACTCTTTTCACATGTT


miRNA
-----

miRNA output are a combination of CircNET and CSCD - for this reason, not all TargetScan or miRanda predictions have complete data available (Energy, Site_type etc.)

.. code-block:: text

    circRNA                     miRNA                  MSA_start  MSA_end  Site_type  Score   Energy  Algorithm
    chr15:100048856-100054054   miR-1208               30.0       37.0     8mer-1a    144.0   -13.62  miRanda,TargetScan
    chr15:100048856-100054054   miR-1237-5p/4488-5p                                   140.0   -20.01  miRanda
    chr15:100048856-100054054   miR-1289                                              141.0   -22.95  miRanda
    chr15:100048856-100054054   miR-141-5p             39.0       45.0     7mer-m8    178.0   -26.76  miRanda,TargetScan

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - circRNA
     - The circRNA coordinates in the hg38 reference genome.
   * - miRNA
     - The miRNA ID.
   * - MSA_start
     - Start coordinate of the miRNA binding site on the circRNA.
   * - MSA_end
     - End coordinate of the miRNA binding site on the circRNA.
   * - Site_type
     - Type of miRNA binding site (8mer-1a, 7mer-m8, 7mer-1a, 6mer).
   * - Score
     - Score of the miRNA binding site.
   * - Energy
     - Minimum free energy of the miRNA binding site.
   * - Algorithm
     - Algorithm(s) used to predict the miRNA binding site (miRanda, PITA, TargetScan).

RBP
---

.. code-block:: text

    circRNA                 Chr     Start       End         RBP_symbol  Information        RBP_numbers  source
    chr18:12493071-12546904 chr18   12503834    12546883    AGO2        SBDH105-59581      24           common
    chr18:12493071-12546904 chr18   12498779    12546872    FBL         SBDH120-145469     10           common
    chr18:12493071-12546904 chr18   12493071    12546757    U2AF2       SBDH19-189369      100          common

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Column Name
     - Description
   * - circRNA
     - The circRNA coordinates in the hg38 reference genome.
   * - Chr
     - Chromosome of the RBP binding site.
   * - Start
     - Start coordinate of the RBP binding site.
   * - End
     - End coordinate of the RBP binding site.
   * - RBP_symbol
     - Symbol of the RNA-binding protein (RBP).
   * - Information
     - Additional information about the RBP binding site.
   * - RBP_numbers
     - Number of RBPs predicted to bind to the circRNA.
   * - source
     - Source of the RBP prediction (cancer, common, normal).