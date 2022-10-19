
# Circall: A fast and accurate methodology for discovery of circular RNAs from paired-end RNA-sequencing data

###########################################################

[![Docker Image CI](https://github.com/datngu/circall/actions/workflows/docker-image.yml/badge.svg)](https://github.com/datngu/circall/actions/workflows/docker-image.yml)


__IMPORTANT NOTE: I tested with 45GB of RAM and 32 CPU cores and recommend a minimum of 45 GB of RAM to obtain stable results for Circall.__

## Update news

### 11 July 2022: version 1.0.1
- Some changes in C++ binary (rebuild the software).

### 06 July 2022: version 1.0.0

- Re-organization of the software package => user-friendly version with Docker.
- Speed up file reading in R with data.table.
- Fix bugs in hg38 annotation.
- Default option "-td" changed to FALSE.
- Output format with an additional column: junction_FPM - junction fragment per million .i.e, junction_fragment_count normalized by million fragment unit.

__NOTE:__ 
- Tutorial has been changed for Docker running tutorial in this version.
- An old version can be found in the sub-directory Circall_v0.1.0.

### 19 June 2020: version 0.1.0

__Instruction how to build the software from source can be found at:__ [https://www.meb.ki.se/sites/biostatwiki/circall/](https://www.meb.ki.se/sites/biostatwiki/circall/)

First submission





## 1. Introduction
Circall is a novel method for fast and accurate discovery of circular RNAs from paired-end RNA-sequencing data. The method controls false positives by two-dimensional local false discovery method and employs quasi-mapping for fast and accurate alignments. The details of Circall are described in its manuscript. In this page, we present the Circall tool and how to use it.

### Software requirements:

- Docker

### Annotation reference

Circall requires:

1) a fasta file of transcript sequences and a gtf file of transcript annotation: can be downloaded from public repositories such as Ensembl (ensembl.org)
2) a genome file of transcript sequences and a gtf file of transcript annotation: can be downloaded from public repositories such as Ensembl (ensembl.org)
3) a RData file of supporting annotation: A description of how to create the RData file for new annotation versions or species is available in the following section.
4) a back-splicing-junction (BSJ) reference sqeuences: A description of how to create the BSJ reference sqeuences for new annotation versions or species is available in the following section.

Current Circall version was tested on the human genome (hg18 and hg38). Detail tutorials are in section 3 (hg19), and section 4 (hg38).


## 2. Download and installation

The esiest way to install and run Circall is using our pre-buit Docker immage: 

```sh
# pull circall docker image:
docker pull ndatth/circall:v1.0.1

# test your docker
docker run --rm ndatth/circall:v1.0.1 Circall.sh
```

Expected output is:

```text
--------------------------------------------------------------------------
                     _
            _____   (_)   _____  _____   __        __      __
           / ___/  / /   / ___/ / ___/  /  \      / /     / /
          / /__   / /   / /    / /__   / __ \    / /__   / /__
          \___/  /_/   /_/     \___/  /_/  \_\  /_/__/  /_/__/

---------------------------------------------------------------------------
                     Checking arguments ....
---------------------------------------------------------------------------


Usage:
./Circall.sh
    -genome [genome in fasta format]
    -gtfSqlite [genome annotation in Sqlite format]
    -txFasta [transcripts (cDNA) in fasta format]
    -txIdx [quasi-index of txFasta]
    -bsjIdx [quasi-index of BSJ reference fasta file]
    -read1 [read1 fastq.gz file]
    -read2 [read2 fastq.gz file]
    -thread [number of thread]
    -tag [tag name]
    -td [TRUE/FALSE value, defaut is FALSE]
    -c[TRUE/FALSE value, defaut is TRUE]
    -o [output folder]

```


__If you want to build Circall from sources: Please prefer to our wiki webpage [https://www.meb.ki.se/sites/biostatwiki/circall](https://www.meb.ki.se/sites/biostatwiki/circall/#user-content-2-download-and-installation) for detail instruction__



## 3. Tutorial for human genome hg19 (GRCh37), annotation version 75.

### 3.1 Prepare BSJ reference database and Sqlite annotation file

#### Download genome fasta, transcript fasta and gtf annotation files.


```sh
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh37.75.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
```

#### Create sqlite

You can generate the sqlite annotation database was generated and able to download from [Homo_sapiens.GRCh37.75.sqlite](https://github.com/datngu/Circall/releases/download/v0.1.0/Homo_sapiens.GRCh37.75.sqlite). This file was generated by the following command:


```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 createSqlite.R \
        data/Homo_sapiens.GRCh37.75.gtf \
        data/Homo_sapiens.GRCh37.75.sqlite
```

#### Create BSJ reference database

The BSJ reference database for Homo_sapiens.GRCh37.75 was generated and able to download from [Homo_sapiens.GRCh37.75_BSJ_sequences.fa](https://github.com/datngu/Circall/releases/download/v0.1.0/Homo_sapiens.GRCh37.75_BSJ_sequences.fa.gz). This file was generated by the following command:


```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 buildBSJdb.R \
        gtfSqlite=data/Homo_sapiens.GRCh37.75.sqlite \
        genomeFastaFile=data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        bsjDist=250 maxReadLen=150 \
        output=data/Homo_sapiens.GRCh37.75_BSJ_sequences.fa
```

### 3.2 Indexing transcriptome and BSJ reference database

#### Index transcriptome

```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
# please uncompress your fa file (with gunzip) before indexing 
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh37.75.cdna.all.fa \
        -o data/IndexTranscriptome
```

#### Index BSJ reference database


```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
# please uncompress your fa file (with gunzip) before indexing 
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh37.75_BSJ_sequences.fa \
        -o data/IndexBSJ
```

Now, you are ready to run Circall.

### 3.3 Run Circall pipeline

Suppose sample_01_1.fasta and sample_01_2.fasta are the input fastq files. For convenience, we prepared a toy example to test the pipeline, which can be downloaded here:

```sh
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_1.fasta.gz
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_2.fasta.gz
```
Circall can be run in one command wrapped in a bash script:

```sh
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 Circall.sh \
    -genome data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    -gtfSqlite data/Homo_sapiens.GRCh37.75.sqlite \
    -txFasta data/Homo_sapiens.GRCh37.75.cdna.all.fa \
    -txIdx data/IndexTranscriptome \
    -bsjIdx data/IndexBSJ \
    -dep Circall/Data/Circall_depdata_human.RData \
    -read1 data/sample_01_1.fasta.gz \
    -read2 data/sample_01_2.fasta.gz \
    -p 4 \
    -tag testing_sample \
    -c FALSE \
    -td FASLE \
    -o data/Testing_HG19_out

```

NOTED: sample depdata of Hela, HEK293, HS68 (cell-lines) are located in Docker immage, path: `Circall/Data/Circall_depdata_human.RData`




## 4. Tutorial for human genome hg38 (GRCh38), annotation version 106.

### 4.1 Prepare BSJ reference database and Sqlite annotation file

#### Download genome fasta, transcript fasta and gtf annotation files.


```sh
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz -O Homo_sapiens.GRCh38.106.gtf.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O Homo_sapiens.GRCh38.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

#### Create sqlite

You can generate the sqlite annotation database was generated and able to download from [Homo_sapiens.GRCh38.106.sqlite](https://github.com/datngu/Circall/releases/download/v1.0.0/Homo_sapiens.GRCh38.106.sqlite). This file was generated by the following command:


```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 createSqlite.R \
        data/Homo_sapiens.GRCh38.106.gtf \
        data/Homo_sapiens.GRCh38.106.sqlite
```

#### Create BSJ reference database

The BSJ reference database for Homo_sapiens.GRCh38.106 was generated and able to download from [Homo_sapiens.GRCh38.106_BSJ_sequences.fa](https://github.com/datngu/Circall/releases/download/v1.0.0/Homo_sapiens.GRCh38.106_BSJ_sequences.fa.gz). This file was generated by the following command:


```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 buildBSJdb.R \
        gtfSqlite=data/Homo_sapiens.GRCh38.106.sqlite \
        genomeFastaFile=data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        bsjDist=250 maxReadLen=150 \
        output=data/Homo_sapiens.GRCh38.106_BSJ_sequences.fa
```

### 4.2 Indexing transcriptome and BSJ reference database

#### Index transcriptome

```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
# please uncompress your fa file (with gunzip) before indexing 
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh38.cdna.all.fa \
        -o data/IndexTranscriptome
```
#### Index BSJ reference database


```sh
# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
# please uncompress your fa file (with gunzip) before indexing 
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh38.106_BSJ_sequences.fa \
        -o data/IndexBSJ
```

Now, you are ready to run Circall.

### 4.3 Run Circall pipeline

Suppose sample_01_1.fasta and sample_01_2.fasta are the input fastq files. For convenience, we prepared a toy example to test the pipeline, which can be downloaded here:

```sh
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_1.fasta.gz
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_2.fasta.gz
```
Circall can be run in one command wrapped in a bash script:

```sh
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 Circall.sh \
    -genome data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -gtfSqlite data/Homo_sapiens.GRCh38.106.sqlite \
    -txFasta data/Homo_sapiens.GRCh38.cdna.all.fa \
    -txIdx data/IndexTranscriptome \
    -bsjIdx data/IndexBSJ \
    -dep Circall/Data/Circall_depdata_human.RData \
    -read1 data/sample_01_1.fasta.gz \
    -read2 data/sample_01_2.fasta.gz \
    -p 4 \
    -tag testing_sample \
    -c FALSE \
    -td FASLE \
    -o data/Testing_HG38_out

```

NOTED: sample depdata of Hela, HEK293, HS68 (cell-lines) are located in Docker immage, path: `Circall/Data/Circall_depdata_human.RData`




## 5. Input, parameters, and output

#### Annotation data

| Parameters      | Description |
| ----------- | ----------- |
| gtfSqlite | genome in fasta format|
| Paragraph | genome annotation in Sqlite format |
| txFasta | transcripts (cDNA) in fasta format |
| txIdx | quasi-index of txFasta|
| bsjIdx | quasi-index of BSJ reference fasta file |


#### Input data

| Parameters      | Description |
| ----------- | ----------- |
| read1 | input read1: should be in gz format|
| read2 | input read1: should be in gz format |



#### Other parameters

| Parameters      | Description |
| ----------- | ----------- |
| dep | data contain depleted circRNAs: to specify the null data (depleted circRNA) for the two-dimensional local false discovery rate method. For convenience, we collect the null data from three human cell lines datasets Hela, Hs68, and Hek293 and provided in the tool: Circall_v0.1.0_linux_x86-64/Data/Circall_depdata_human.RData |
| p | the number of threads: Default is 4 |
| tag | tag name of results: Default is “Sample” |
| td | generation of tandem sequences: TRUE/FALSE value, default is FASLE |
| c | clean intermediate data: TRUE/FALSE value, default is TRUE |
| o | output folder: Default is the current directory |


#### Output

The main output of Circall is provided in \**_Circall_final.txt.* In this file, each row indicates one circular RNA, and the information of one circular RNA is presented in 8 columns:

| Name      | Description |
| ----------- | ----------- |
| chr | chromosome |
| start | start position |
| end | end position |
| geneID | gene name that the circRNA belongs to |
| circID | the ID of circRNA in the format “chr__start__end” |
| junction_fragment_count | the number of fragment counts supporting the back-splicing-junction (BSJ) |
| junction_FPM | junction fragment per million .i.e, junction_fragment_count normalized by million fragment unit |
| median_circlen | the median length of the circular RNA |
| fdr | the false discovery rate computed from the two-dimensional local false discovery method |




## 6. Circall simulator

### Introduction

Circall simulator is a tool integrated in Circall to generate RNA-seq data of both circRNA and tandem RNA. The source codes are provided in `R/Circall_simulator.R` of the Circall package. The main function of the simulator is `Circall_simulator()` which is able to be run in R console. This function requires the following parameters:

### Parameter setting:

| Name      | Description |
| ----------- | ----------- |
| circInfo | a data frame that contains 6 columns which are: Chr, start_EXONSTART, end_EXONEND, GENEID, cCount and FPKM. Chr is chromosome name with formated as 1:22, X, Y, Mt. start_EXONSTART is starting position starting exon of circRNA, end_EXONEND is ending position ending exon of circRNA, GENEID is gene ID contains circRNA (used to get gene model), cCount are number of read pair want to generate for the target circRNA and FPKM are Fragments Per Kilobase of transcript per Million of target circRNAs. This is used to simulate circular RNAs|
| tandemInfo | a data frame similar to circInfo to simulate tandem RNAs. tandemInfo=NULL (the default value) to not simulate tandem RNAs|
| error_rate | sequencing error rate, the default value is 0.005|
| set.seed | set seed for reproducibility, the default value is 2018|
| gtfSqlite | path to your annotation file, Sqlite formated (generated by GenomicFeatures)|
| genomeFastaFile | path to your genome fasta file|
| txFastaFile | path to your transcript fasta file (cDNA)|
| out_name | prefix output folders, the default value is “Circall_simuation”|
| out_dir | the directory contains output, the default value is the current directory|
| lib_size | expected library size used when useFPKM=TRUE, the default value is NULL|
| useFPKM | boolean value to use FPKM or not, the default value is FALSE. When this useFPKM=TRUE, users need to set value for lib_size, and the simulator will use the abundance in column FPKM of circInfo/tandemInfo for simulation|

### A toy example for using Circall simulator

For an illustration of using Circall simulator, we provide in this section a toy example. Suppose your current working directory contains the installed Circall and the annotation data. First, we need to load the functions of the simulator into your R console:

```R
source("/Circall_v1.0.0/R/Circall_simulator.R")
```
Then we create objects __circInfo__ and __tandemInfo__ containing the information of CircRNAs and tandem RNAs

```R
Chr = c(7,7,3,5,17,4,7,1,3,1,17,12,14,10,18,17,5,20,16,17)

start_EXONSTART = c(131113792,99795401,172363413,179296769,36918664,151509200,2188787,51906019,57832924,225239153,76187051,111923075,104490906,101556854,196637,21075331,74981032,60712420,56903641,80730328)

end_EXONEND = c(131128461,99796580,172365904,179315312,36918758,151509336,2270359,51913807,57882659,225528403,76201599,111924628,104493276,101572901,199316,21087123,74998635,60716000,56904648,80772810)

GENEID = c("ENSG00000128585","ENSG00000066923","ENSG00000144959","ENSG00000197226","ENSG00000108294","ENSG00000198589","ENSG00000002822","ENSG00000085832","ENSG00000163681","ENSG00000185842","ENSG00000183077","ENSG00000204842","ENSG00000156414","ENSG00000023839","ENSG00000101557","ENSG00000109016","ENSG00000152359","ENSG00000101182","ENSG00000070915","ENSG00000141556")

set.seed(2021)
cCount = sample(2:2000,20)
FPKM = rep(0,20)

BSJ_info = data.frame(Chr = Chr, start_EXONSTART = start_EXONSTART, end_EXONEND = end_EXONEND, GENEID = GENEID, cCount = cCount, FPKM = FPKM)

circSet=c(1:15)
circInfo = BSJ_info[circSet,]
tandemInfo = BSJ_info[-circSet,]
```

Finally, we run the simulator:

```R
simulation = Circall_simulator(circInfo = circInfo, tandemInfo = tandemInfo, useFPKM=FALSE, out_name = "Tutorial", gtfSqlite = "Homo_sapiens.GRCh37.75.sqlite", genomeFastaFile = "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa", txFastaFile = "Homo_sapiens.GRCh37.75.cdna.all.fa", out_dir= "./simulation_test")

```

You can find in the “./simulation_test” that contains the outputs including:

- simulation_setting: setting information of simulation of both circRNAs and tandem RNAs.
- circRNA_data: RNA seq data of CircRNAs
- tandem_data: RNA seq data of tandem RNA
- fasta sequences of tandem RNAs
- fasta sequences of circular RNAs


## 7. License

Circall uses GNU General Public License GPL-3.

## 8. Reference and citation

Nguyen, Dat Thanh, Quang Thinh Trac, Thi-Hau Nguyen, Ha-Nam Nguyen, Nir Ohad, Yudi Pawitan, and Trung Nghia Vu. 2021. “Circall: Fast and Accurate Methodology for Discovery of Circular RNAs from Paired-End RNA-Sequencing Data.” BMC Bioinformatics 22 (1): 495. https://doi.org/10.1186/s12859-021-04418-8.
