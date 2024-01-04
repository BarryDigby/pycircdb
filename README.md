# pycircdb

## Reference files used

You must enforce these when creating each database pair-key values. 

### hg19

* `Cerina` (contains circAtlas, circBase)
* `Circbase`
* `CircBank`
* `Arraystar`

```console
> circbase
chr1    1158623 1159348 -       hsa_circ_0000002

> circbank
hsa_circSDF4_003        hsa_circ_0000002        chr1:1158623-1159348

> cerina
chr1:1158623-1159348    hsa-SDF4_000001 hsa_circ_0000002
```

circbase, circbank and cerina (column 1 not CIRI2) are all 0-based coutning.

### hg38

* `circatlas`
* `cscd`
* `circnet`

> CSCD is 0-based, circAtlas is 1-based

not appropriate to convert circBank identifiers using chain - the underlying gene may change when converting references. 

```console
> cscd
chr5:136152674|136154163    hsa_circ_0074050

> circatlas
chr5    136152675   136154163  +   hsa-SMAD5_0007
```

make decision to make all 0-based or all 1-based.

# miRNA targets 

again figure out the reference genomes used here. 

***

unique circbase ids in cscd: 86652
unique circbase id in hg19_IDs that are in cscd: 85845

so all that is needed is to merge these missing position + circbase IDS to the hg19_lift_IDs.txt file...

Then you have a proper liftover of hg19 and hg38. 


## Databases

### circ2disease

Missing circRNA-miRNA information for download.

* `all_disease.txt`: List of diseases in database
```console
abdominal aortic aneurysm
acoustic neuroma
acquired immunodeficiency syndrome
```
* `circ2disease_association.txt`: Main database file containing `circBase name`, `Synonymns`, `Chromosome`, `Start-End`, `Strand`, `Genome version`, `Downstream gene`, `Sequence`, `circRNA type`, `Disease`, `Method(s) found`, `Sample(s) found`, `Cellular location`, `Expression pattern`, `miRNA validated`, `miRNA target validated`, `RBP validated`, `RBP target validated`, `Functional description`, `PubMed ID`, `Year`, `Journal`, `Title`, `Note`. Tab separated.

only 273 entries. Majority are hg19 or N/A. Only 4 GRCh38 entries.

```console
circLARP4       hsa_circ_0003222; hsa_circRNA_101057    12      50848096-50855130       +       hg19    LATS1   GCTCCCTTTCCCAATGGTAGTTTTGTGAATGGCTTTAATTCGCCAGGATCTTATAAAACAAATGCTGCTGCTATGAATATGGGTCGACCATTCCAAAAAAATCGTGTGAAGCCTCAGTTTAGGTCATCTGGTGGTTCAGAACACTCAACAGAGGGCTCTGTATCCTTGGGGGATGGACAGTTGAACAGATATAGTTCAAGAAACTTTCCAGCTGAACGGCATAACCCCACAGTAACTGGGCATCAGGAGCAAACTTACCTTCAGAAGGAGACTTCCACTTTGCAGGTGGAACAGAATGGGGACTATGGTAGGGGCAG  exonic   gastric cancer  Microarray; in vitro knockdown; in vitro overexpression; dual luciferase reporter assay;  RIP; FISH; etc.       GES-1 and GC cell lines (SGC-7901, MKN-45, MKN-28, HGC-27, MGC-803, AGS, BGC-823)       cytoplasm       down-regulated  miR-424-5p    LATS1   N/A     N/A     Increased miR-424 expression or decreased LATS1 expression was associated with pathological stage and unfavorable prognosis of GC patients. Ectopic expression of miR-424 promoted proliferation and invasion of GC cells by targeting LATS1 gene. Furthermore, circLARP4 was mainly localized in the cytoplasm and inhibited biological behaviors of GC cells by sponging miR-424. The expression of circLARP4 was downregulated in GC tissues and represented an independent prognostic factor for overall survival of GC patients.   28893265        2017    Molecular Cancer       Circular RNA_LARP4 inhibits cell proliferation and invasion of gastric cancer by sponging miR-424-5p and regulating LATS1 expression     N/A
```

* `circ_disease.txt`: List of diseases
```console
acute myeloid leukemia
Alzheimer's disease
angiogenesis
atherosclerosis
basal cell carcinoma
```
* `circ_group.txt`: circRNA ID and corresponding synonyms. Tab separated.

Only 237 entries

```console
circHIPK3       BCRC-2  circR-284       hsa_circRNA_100782      hsa_circ_000016 hsa_circ_0000284
cAFF1   hsa_circ_0001423
cDENND4C        circDENND4C     hsa_circRNA_005684      hsa_circ_0005684
```
* `circRNA_RBP_predict.txt`: Predicted RBP(#number binding sites) - `circBase name`, `(first) synonym`, `RBP`.

Only 238 entries

```console
CircRNA_Name    circBase        RBP(#Number of binding sites)
circHIPK3       hsa_circ_0000284        EIF4A3(15); AGO2(12); FMRP(7); PTB(5); FUS(4); LIN28A(3); IGF2BP3(3); HuR(3); IGF2BP2(2); MOV10(1); IGF2BP1(1); CAPRIN1(1); AGO1(1)
cAFF1   hsa_circ_0001423        EIF4A3(10); AGO2(7); PTB(3); IGF2BP3(3); IGF2BP1(3); IGF2BP2(2); HuR(2); FMRP(2); FUS(1); DGCR8(1); C17ORF85(1)
```
* `HMDD_and_OncomiRDB_and_dbDEMC_filter.txt`: miRNAs validated by one of three databases. `miRNA ID`, `Disease`, `Description (expression, recurrence etc)`, `PMID`, `Database derived from`.
```console
hsa-miR-206     muscular dystrophies    miR-206: miR-206 high expression for muscle regeneration and maturation 18827405        HMDD v2.0
hsa-miR-182-5p  ovarian cancer  up      24512620        dbDEMC v2.0
```
* `miRecords_and_miRTarBase_filter.txt`: Expreimentally validated miRNAs. `miRNA ID`, `Target gene`, `Validation assay`, `PMID`.
```console
hsa-miR-20a-5p  HIF1A   Luciferase reporter assay; Western blot; Northern blot; qRT-PCR 18632605
hsa-miR-20a-5p  HIF1A   Luciferase reporter assay; qRT-PCR; Western blot        23911400
hsa-miR-21-5p   RASGRP1 qRT-PCR; Luciferase reporter assay; Western blot        20483747
```

### circatlas

* `human_bed_v3.0.txt`:

```console
Chro    Start   End     Strand  circAltas_ID
KI270742.1      1326    34982   -       hsa-intergenic_003508
chr1    805799  810170  -       hsa-RP11-206L10_0001
chr1    808574  810170  -       hsa-RP11-206L10_0002
chr1    926703  927095  -       hsa-SAMD11_0013
```

### circbank

* `circBank_circrna_annotation.txt`

```console
circBankID      circbaseID      position        strand  splicedSeqLength        annotation      bestTranscript  geneSymbol
hsa_circA1CF_001        hsa_circ_0018410        chr10:52575765-52623840 -       1470    ANNOTATED, CDS, coding, INTERNAL, OVCODE, OVERLAPTX, OVEXON, UTR5       NM_001198819    A1CF
hsa_circA1CF_002        hsa_circ_0018409        chr10:52559168-52573822 -       7964    ANNOTATED, CDS, coding, OVCODE, OVERLAPTX, OVEXON, UTR3 NM_001198819    A1CF
```

* `circbank_miRNA_all_v1.txt`

```console
miR_ID  circbase_ID     Tot Score
hsa-miR-6087    hsa_circ_0002172        11645.00
hsa-miR-6087    hsa_circ_0089776        11645.00
```

* `circRNA m6a0528.txt`

```console
chr     start   end     strand  name    Read counts in Eluate   Read counts in supernatant      m6A levels      sample
chr2    100046301       100052403       -       chr2junc7192_R1_chr2junc92_A    2       0.01    0.99227391       hESCs
chr2    100055062       100081447       -       chr2junc75_R1_chr2junc113_A     2       0.01    0.99227391       hESCs
```

* `coding potential/CPATtriResult.txt`

```console
circBankID      circRNA_size    ORF_size        Fickett_score   Hexamer_score   coding_prob
hsa_circ_chr1_00102     321     111     0.5331  0.0145947351656 0.00921202797119123
hsa_circSDF4_003        753     222     1.2124  0.580393899499  0.910132807549632
hsa_circGNB1_039        681     231     0.7451  0.0874645903186 0.101360579092505
```

* `coding potential/IRESpoition.txt`

```console
id      start   end     score
hsa_circ_0000003        2051    2224    0.9473 
hsa_circ_0000003        4451    4624    0.9570 
hsa_circ_0000003        7101    7274    0.9745 
hsa_circ_0000003        8851    9024    0.9647 
```

* `hg19_circbase_seq.fa`

```console
>hsa_circ_0000001
ATGGGGTTGGGTCAGCCGTGCGGTCAGGTCAGGTCGGCCATGAGGTCAGGTGGGGTCGGCCATGAAGGTGGTGGGGGTCA
TGAGGTCACAAGGGGGTCGGCCATGTG
>hsa_circ_0000002
GGTGGATGTGAACACTGACCGGAAGATCAGTGCCAAGGAGATGCAGCGCTGGATCATGGAGAAGACGGCCGAGCACTTCC
AGGAGGCCATGGAGGAGAGCAAGACACACTTCCGCGCCGTGGACCCTGACGGGGACGGTCACGTGTCTTGGGACGAGTAT
AAGGTGAAGTTTTTGGCGAGTAAAGGCCATAGCGAGAAGGAGGTTGCCGACGCCATCAGGCTCAACGAGGAACTCAAAGT
GGATGAGGAAA
```

## circnet

* `{disease}.csv`

```console
,CircID,Strand,CircType,HostGene,Algorithm,Sequence
0,chr7:139715932|139717015,-,exon,ENSG00000064393.15;,find_circ;CIRI2;DCC;,GTATGGCCTCACATGTGCAAGTTTTCTCCCCTCACACCCTTCAATCAAGTGCCTTCTGTAGTGTGAAGAAACTGAAAATAGAGCCGAGTTCCAACTGGGACATGACTGGGTACGGCTCCCACAGCAAAGTGTATAGCCAGAGCAAGAACATCCCCCTGTCGCAGCCAGCCACCACAACCGTCAGCACCTCCTTGCCGGTCCCAAACCCAAGCCTACCTTACGAGCAGACCATCGTCTTCCCAGGAAGCACCGGGCACATCGTGGTCACCTCAGCAAGCAGCACTTCTGTCACCGGGCAAGTCCTCGGCGGACCACACAACCTAATGCGTCGAAGCACTGTGAGCCTCCTTGATACCTACCAAAAATGTGGACTCAAGCGTAAGAGCGAGGAGATCGAGAACACAAGCAGCGTGCAGATCATCGAGGAGCATCCACCCATGATTCAGAATAATGCAAGCGGGGCCACTGTCGCCACTGCCACCACGTCTACTGCCACCTCCAAAAACAGCGGCTCCAACAGCGAGGGCGACTATCAGCTGGTGCAGCATGAGGTGCTGTGCTCCATGACCAACACCTACGAGGTCTTAGAGTTCTTGGGCCGAGGGACGTTTGGGCAAGTGGTCAAGTGCTGGAAACGGGGCACCAATGAGATCGTAGCCATCAAGATCCTGAAGAACCACCCATCCTATGCCCGACAAGGTCAGATTGAAGTGAGCATCCTGGCCCGGTTGAGCACGGAGAGTGCCGATGACTATAACTTCGTCCGGGCCTACGAATGCTTCCAGCACAAGAACCACACGTGCTTGGTCTTCGAGATGTTGGAGCAGAACCTCTATGACTTTCTGAAGCAAAACAAGTTTAGCCCCTTGCCCCTCAAATACATTCGCCCAGTTCTCCAGCAGGTAGCCACAGCCCTGATGAAACTCAAAAGCCTAGGTCTTATCCACGCTGACCTCAAACCAGAAAACATCATGCTGGTGGATCCATCTAGACAACCATACAGAGTCAAGGTCATCGACTTTGGTTCAGCCAGCCACGTCTCCAAGGCTGTGTGCTCCACCTACTTGCAGTCCAGATATTACAG
```

* `circRNA_miRNA interaction.csv`

```console
CircID,mir1,PITA,miRanda,targetScan
chr3:195897004|195897224,hsa-let-7b-5p,0.0,0.0,1.0
chr5:151673190|151674644,hsa-let-7b-5p,0.0,0.0,1.0
chr12:118033598|118033703,hsa-let-7b-5p,0.0,0.0,1.0
```

* `Samples.csv`

```console
circ_id,sample_source,sample_type,sample_id,read
chr7:139715932|139717015,mioncocirc,prostate,VCaP-capt-SI_8045-C4D5HACXX,220
chr7:139715932|139717015,mioncocirc,prostate,VCaP-capt-SI_8046-C4D5HACXX,235
chr7:139715932|139717015,mioncocirc,prostate,VCaP-capt-SI_8047-C4D5HACXX,231
```

## CSCD 

Contains circ annotation, miRNA MRE sites, RBP sites for cancer, normal and common.

* `hg38_normal_circrna_circ_chrY.txt`

```console
circRNA circbase_id     type    sample_source   position        host_gene       reads_counts    algorithm       chain   subcellular_location    gene_type       average_read_counts     average_log2SRPTM       log2SRPTM       alternative_splicing    ratios  works_cited     mioncocirc      number_of_algorithm     sequence
chrY:2797471|2800763    ---     normal  Laryngeal_mucosa_22     Intergenic      ---     Laryngeal_mucosa_22,find_circ,1 find_circ       +       ---     ---     1.0     8.297151336774899       Laryngeal_mucosa_22,find_circ,8.297151336774899 ---     ---     ---     ---     1       Laryngeal_mucosa_22,TGCTGGGAATACAGATATGAGCCACCACGCCCGACCTTCTAAAGTAAAAAGAGTTACAGTGAGCTAAGATTAATTTGTTGTTGATTAACAATTACATTTAATGACTGGGCCCAGTGGCTAATGCCTGTAATCCCAGCACTTTAGGAGACCAAGGCAGGCAGATAACTTGAGGTCAGGAGTTGAAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAGCACAAAAATTAGCTGGGCATGTTGGCTCACGCCTGTAATCCCAGCTACTTGGGAGACTGAGGCAGGAGAATCACTTGAACTTGGGAGGCAGAGGTTACAGTTAGACAAGATGGCATCACTGCAATCCAGCCTGGGTGACAGAGCGAGACTCTGTCTCAAAATAAAAAAATAAACATAAATTTAGTGTAGCCTAAGTGTTTGTGAAGCATACAGTATTACACAGTAGGGTCCTGGGCCTTCACCTTCACTCACTACTCACTCACTGGCTCACCCAGAGCAACGTCCAGCCCGCAAGTGCTAATCATGGCAAGTATCCTACCTGTATAGCTGTATTACTTTTTTTCATGTTTTCTACTATATTTTTTACTGTATATTTTCCATGTTTACACACAAATACTCACCATTGGGTTACAATTGCCAACTGGATTTCTTGCAGTAACTTGCTCTACAGGTTTGTAGCCTAGGAGTGCTAGGCTCTACCATACAGCCTAGGTGTGTAATAGGCTACACCCCCTAGGTTTCTGTAAATACACTTTATCATGTGCTTACAATGAGGAAATCACCTAATGATGCATTTCTCAGTAGGTATCCCCATTATTGAGTGATTGATGCAGGACTATAGCTTAAATGTGAATTCATTTGAAGTAGTTTGCATGATTAGCAGTTTGATATCTGGATTGTCTTTTATTCTTGTTAACACACACACAAATACACACACACACACACACACACACACTCCTGTTAGATAAAAACCTTGTCTCTGGCTGAGTGCGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGGGGCCAAGGCAGGGAGATCACAAGGGCAGAAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCTGGGCGTGTGTGCACACCTGTAATCCCAGCTATTTGGGAGGCTGAGGCAGAATTGCTTCAACCTGGGAGGCAGAGGTTACAGTGAACCAAAATTGCACCACTGCATTCCAGTCTGGGTAACAGAGCAAGACTCCATCTCAAAAACACAAAACAAAACAAAACAAAAAACTGTCTCTGTCTATATTTTGTCTGGGCTCACCTTTTCTTATTTGAGATGCTTATAAATTCATTGTGTGATTTTGTTTTGCTCTGTTGAAAATGTGAAGATGAATCTCTGAATAAAAACGTTTTATTTGGGAAGTAAGACTTGCAATCTGGTGCCTACATGCAGACAGTGTGGTCATCAGTATGTCTGATGAACAAACACAAGCTTTAAGATTTTATAAAAGAAAAAGAGAAATGGTGTATTGTTCCTTGAGAAAGTTCATTGGCACTATTCAGGTTTTGGAGAGCTGGCAAGCAAGCTGATGGTGGTGCATAAAACTAGTCTTAGATTTTTAGCAAGTTGTTTCAGTAGCCAGGCTAGCTGAGAGTTACATTCTTGGAACAATCTTTTGTTGCTGTTGTTGTTTTGAGAGAGAGTCTAGCTCTGTCCCCAGGCTGGAGTGCAGTGGTGCAATTTTGGCTCACTGCAACCTCTGCCTCCCAGGTTCAAGCAATTCTGCCTCAGCCTCTCAAGTAGCTGGGATTACACGTGTGTGCATCCATGCCCAGCTAATTTTTGCATTTTTAGTAGAGACAGGGTTTCCCCATGTTGGCCAGGCTGGTCTCAAACTCCTAGCCTCAAGTGATCCGCCAGCCTTGGCTTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCACAGCCAGTCTGGAACAATGTTATGTAACCTGTTTCCTCACCCATTATGGCCCCACAAACCCCAGCCTTTCAACTCTGTTTCAGTTGGCTATCACGAGTGACCCAACTTGTATAATCAACTTTTACGTTTTTCCATTTTGATAAAGACCTTTCTCTGAAAGCATTGCTGATCAACCATTTTGAAGTTAGGCTTAATTGTCCATTAATATCTGGATGGACCTGTCCTGGTTGCTTTCAGTCATGTTGAGGGGAACATTATGGAGGTCACATTTGAGGCCAGAAGGGAACATTTTCTCAAGCAAATCTGTCTGGAGTCCAGCATCAAGTTCTGTCTTATTAGTCGATAAGCATCAGCAACTGTCTCAAAGTGTTGCATTAACTTTATCTTTTTCTTTTGTTGAAATGGAGTTTCACTCTTGTTGCCCAGGCTGGAGTGTAATGGCACAATCTTGGCTCACTGAAACCTGCACCTCCCAGGTTCAAGTGATTCTCCTGCTTCACCCTCCCAAGTAGTTGGGATTACAGGCATGATCCACTATGCCCAGCTAATTTTTTGTATTTAGTAAAAATGGTTTCACCATGTTTGTCAAGCTGGTCTTGAATTCACTCCTGACCTCAGGTGATCCACCCTCCTCAGCCTCCCAAAGTGCTGGGATTGCAGGCGTGAGCCACTGCACCCGAGAACTTTATCTTTGTAGGAGAGTTGGCTTTACGAGGATTAGACAGGCAATAAAAACAAAGTTTAAAAAAATATATAAAAATAATTAATAGTAATATTATAAATTTAGTTTGTTTAACGGTTTTAAACTAAGATCCTAAGACTAAGAGCAAAAAGCCAGACAAATAAAAAGCCTGTGGGGTAAATTGGATGAAATCTGTGTTAGTTGTTTAGTAGATCTTATGACTAGATTGACTTAAAGCAAAGAACACCAGCTTTATAAAAAGGGACCCACTAGTATGATTTGAACAAAAAGTTGTGATACTGGTATTTTCAGTTGGCCACTAGGTGGACTAAAAGGATTCCTTATGCCAGGTTTTGTTGAATTACCGATAGCAGTTACTGACTGGAAAATTTTATCCTGCCAAATGAAGAAGAAGGTAGCATTAGCGGGGTGATAGTCTCATTCTGATTTAAAGTCTTGTTCTGATGTCTTAGGAAAAGTTCTGTACAGTGTGGGAAAATGTCAACTTTTTGTCTTGTCTTGAATGTCTCTGGTGATGGCATTGGGTGGTTTGGTAAACTCTTGGTGCAGCTAATGCAAATCTATTTATTTTGTGATTTGCTGCCTCATGTATCATGCAGAAGACTTGTTTTTAAAAACATGT
```

* `hg38_normal_circrna_MRE_chrY.txt`

```console
CircRNA_ID      microRNA_ID     MSA_start       MSA_end Site_type       Score   Energy  Align_start     Align_end       Algorithm       Number_of_algorithm     Numbers
chrY:1880496|1888746    miR-1228-3p     71      76      6mer    ---     ---     ---     ---     TargetScan      1       1
chrY:1880496|1888746    miR-124-3p.1    34      40      7mer-m8 151.00  -17.14  19      41      TargetScan,miRanda      2       1
```

* `hg38_normal_circrna_RBP_chrY.txt`

```console
CircRNA_ID      Chr     Start   End     RBP_symbol      Information     RBP_numbers
chrY:2797471|2800763    chrY    2797613 2797614 UPF1    SBDH155-101630  1
chrY:2797471|2800763    chrY    2799208 2799227 PTBP1   SBDH22-51266    1
chrY:2797471|2800763    chrY    2799493 2799511 ELAVL1  human_RBP_CLIPdb_28753945;GSE50989,GSM1234284   1
chrY:2797471|2800763    chrY    2800019 2800079 AGO1-4  SBDH95-122411   1
```

## exorbase

* `circRNAs_anno.txt`

```cosnole
circID  circBase ID     Genomic position        Strand  Genomic length  Spliced length  Gene symbol     Gene type       Sample type     Tumor frequency (Sample number) Benign frequency (Sample number)        Healthy frequency (Sample number)       Urine frequency (Sample number) CSF frequency (Sample number)   Bile frequency (Sample number)  Tumor mean      Benign mean     Healthy mean    Urine mean      CSF mean        Bile mean       Diff group      Plot
exo_circ_00001  hsa_circ_0017425        chr10:1000677-1001013   +       336     256     GTPBP4  protein coding gene     Urine   0(0)    0.008(1)        0(0)    0.062(1)        0(0)    0(0)    0       0.04    0      5.26     0       0       NA      NA
exo_circ_00002  hsa_circ_0017426        chr10:1000677-1007128   +       6451    455     GTPBP4  protein coding gene     Blood   0.006(4)        0.008(1)        0(0)    0(0)    0(0)    0(0)    3.12    0.04    0      00       0       NA      NA
```

* `{Disease}_circRNA.txt`: Expression matrices of exo_circrnas in samples.

## mioncocirc

* `v0.1.release.txt`

```console
chr     start   end     reads   symbol  sample  release
chr7    139715931       139717015       220     HIPK2   VCaP-capt-SI_8045-C4D5HACXX     v0.0
chr21   36247516        36248568        175     DOPEY2  VCaP-capt-SI_8045-C4D5HACXX     v0.0
chr5    137985256       137988315       156     FAM13B  VCaP-capt-SI_8045-C4D5HACXX     v0.0
chr9    135881632       135883078       154     CAMSAP1 VCaP-capt-SI_8045-C4D5HACXX     v0.0
chr1    117402185       117420649       139     MAN1A2  VCaP-capt-SI_8045-C4D5HACXX     v0.0
```