
# Documents for Circall version 0.0.0

###########################################################

## Update news


#### 30 May 2020: version 0.0.0
- First submission

## 1. Introduction
Circall is a novel method to discover circular RNA from paired-end RNA sequencing data. Circall is characterized by employing quasi-mapping for fast and accurate alignments and the multidimensional local false discovery method for circRNA candidate assessments that improve circRNA detection accuracy. Full details of the method is described in the method publication.

### Software requirements:
Circall is implemented in R and C++. We acknowledge for materials from Sailfish, Rapmap and other tools used in this software.
- A C++-11 compliant compiler version of GCC (g++ >= 4.8.2)
- R packages version 3.6.0 or latter with following installed packages: GenomicFeatures, Biostrings, foreach, and doParallel.

### Annotation reference
Circall requires

1) a fasta file of transcript sequences and a gtf file of transcript annotation: can be downloaded from public repositories such as Ensembl (ensembl.org)
2) a genome file of transcript sequences and a gtf file of transcript annotation: can be downloaded from public repositories such as Ensembl (ensembl.org)
3) a RData file of supporting annotation: A description of how to create the RData file for new annotation versions or species is available in Section X.

Current Circall version was tested on the human genome, transcriptome with ensembl annotation version GRCh37.75. Following files are required:
- Sequences of genome (ensembl website) [GRCh37.75 genome fasta](http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz)
- Sequences of transcripts (ensembl website) [GRCh37.75 cdna fasta](http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz)
- Gtf annotation of transcripts (ensembl website) [GRCh37.75 gtf annotation](http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)

### Versions
The latest version and information of Circall is updated at: https://github.com/datngu/Circall

The older versions can be found here:
- Version 0.0.0: https://github.com/datngu/Circall/releases/tag/v0.0.0

## 2. Download and installation
If you use the binary verion of Circall: 
- Download the lastest binary version from Circall website: [Circall_v0.0.0](https://github.com/datngu/Circall)
```sh
wget https://github.com/datngu/Circall/releases/download/v0.0.0/Circall_v0.0.0_linux_x86-64.tar.gz -O Circall_v0.0.0_linux_x86-64.tar.gz
```
- Uncompress to folder
```sh
tar -xzvf Circall_v0.0.0_linux_x86-64.tar.gz
```
- Move to the *Circall_home* directory and do configuration for Circall
```sh
cd Circall_v0.0.0_linux_x86-64
bash config.sh
cd ..
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/Circall_v0.0.0_linux_x86-64/linux/lib:$LD_LIBRARY_PATH
export PATH=/path/to/Circall_v0.0.0_linux_x86-64/linux/bin:$PATH
```
### Do not forget to replace "/path/to/" by your local path.

or used this command to automaticlly replace your path:

```sh
export LD_LIBRARY_PATH=$PWD/Circall_v0.0.0_linux_x86-64/linux/lib:$LD_LIBRARY_PATH
export PATH=$PWD/Circall_v0.0.0_linux_x86-64/linux/bin:$PATH
```

If you want to build Circall from sources:

write later...


## 3. Prepare BSJ reference database and Sqlite annotation file
### Download genome fasta, transcript fasta and gtf annotation files.
```sh
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh37.75.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
```

### Create sqlite
```sh
Rscript Circall_v0.0.0_linux_x86-64/R/createSqlite.R Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.sqlite
```
### Create BSJ reference database
```sh
Rscript Circall_v0.0.0_linux_x86-64/R/buildBSJdb.R gtfSqlite=Homo_sapiens.GRCh37.75.sqlite genomeFastaFile=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa bsjDist=250 maxReadLen=150 output=Homo_sapiens.GRCh37.75
```
## 4. Indexing transcriptome and BSJ reference database
### Index transcriptome
```sh
Circall_v0.0.0_linux_x86-64/linux/bin/TxIndexer -t Homo_sapiens.GRCh37.75.cdna.all.fa -o IndexTranscriptome
```
### Index BSJ reference database

```sh
Circall_v0.0.0_linux_x86-64/linux/bin/TxIndexer -t Homo_sapiens.GRCh37.75_BSJ_sequences.fa -o IndexBSJ
```
Now, you are ready to run Circall.

## 5. Run Circall pipeline

You can download a prepared test data to test the pipeline:
```sh
wget https://github.com/datngu/Circall/releases/download/v0.0.0/sample_01_1.fasta.gz
wget https://github.com/datngu/Circall/releases/download/v0.0.0/sample_01_2.fasta.gz
```

You can run Circall in one command that is warpped as a bash script:

```sh
bash Circall_v0.0.0_linux_x86-64/Circall.sh -genome Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -gtfSqlite Homo_sapiens.GRCh37.75.sqlite -txFasta Homo_sapiens.GRCh37.75.cdna.all.fa -txIdx IndexTranscriptome -bsjIdx IndexBSJ -dep Circall_v0.0.0_linux_x86-64/Data/Circall_depdata_human.RData -read1 sample_01_1.fasta.gz -read2 sample_01_2.fasta.gz -p 4 -tag testing_sample -c FALSE -o Testing_out
```
### Obligatory parameters are:
- genome -- genome in fasta format
- gtfSqlite -- genome annotation in Sqlite format
- txFasta -- transcripts (cDNA) in fasta format
- txIdx -- quasi-index of txFasta
- bsjIdx -- quasi-index of BSJ reference fasta file
- read1 -- input read1: should be in gz format
- read2 -- input read2: should be in gz format

### Optional parameters are:
- dep -- data contain depleted circRNAs: specify location of deletep data that is used as training data for estimation of fdr2d. We have prepared a Rdata file that contain data dirived from Hela, Hs68, and Hek293 datasets: Circall_v0.0.0_linux_x86-64/Data/Circall_depdata_human.RData
- p -- number of thread: Defaut is 4
- tag -- tag name of results:  Defaut is "Sample"
- td -- generation of tandem sequences: TRUE/FALSE value, defaut is TRUE
- c -- clean intermediate data: TRUE/FALSE value, defaut is TRUE
- o -- output folder: Defaut is current directory




## 6. A practical copy paste example of HEK293 dataset

### Download and install Circall
```sh
wget https://github.com/datngu/Circall/releases/download/v0.0.0/Circall_v0.0.0_linux_x86-64.tar.gz -O Circall_v0.0.0_linux_x86-64.tar.gz
```
- Uncompress to folder
```sh
tar -xzvf Circall_v0.0.0_linux_x86-64.tar.gz
```
- Move to the *Circall_home* directory and do configuration for Circall
```sh
cd Circall_v0.0.0_linux_x86-64
bash config.sh
cd ..
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH

```sh
export LD_LIBRARY_PATH=$PWD/Circall_v0.0.0_linux_x86-64/linux/lib:$LD_LIBRARY_PATH
export PATH=$PWD/Circall_v0.0.0_linux_x86-64/linux/bin:$PATH
```

### Download genome fasta, transcript fasta and BSJ databases and annotation file.

```sh
# genenome
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
# cDNA (transcript)
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh37.75.cdna.all.fa.gz
# pre-built BSJ databases
wget https://github.com/datngu/Circall/releases/download/v0.0.0/Homo_sapiens.GRCh37.75_BSJ_sequences.fa.gz
gunzip Homo_sapiens.GRCh37.75_BSJ_sequences.fa.gz
# pre-genarated Sqlite annotation
wget https://github.com/datngu/Circall/releases/download/v0.0.0/Homo_sapiens.GRCh37.75.sqlite
```
### Index transcriptome
```sh
Circall_v0.0.0_linux_x86-64/linux/bin/TxIndexer -t Homo_sapiens.GRCh37.75.cdna.all.fa -o IndexTranscriptome
```
### Index BSJ reference database

```sh
Circall_v0.0.0_linux_x86-64/linux/bin/TxIndexer -t Homo_sapiens.GRCh37.75_BSJ_sequences.fa -o IndexBSJ
```

### Download HEK293 RNA seq data

```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR347/003/SRR3479243/SRR3479243_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR347/003/SRR3479243/SRR3479243_2.fastq.gz
```

### Run Circall

```sh
bash Circall_v0.0.0_linux_x86-64/Circall.sh -genome Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -gtfSqlite Homo_sapiens.GRCh37.75.sqlite -txFasta Homo_sapiens.GRCh37.75.cdna.all.fa -txIdx IndexTranscriptome -bsjIdx IndexBSJ -dep Circall_v0.0.0_linux_x86-64/Data/Circall_depdata_human.RData -read1 SRR3479243_1.fastq.gz -read2 SRR3479243_2.fastq.gz -p 4 -tag testing_sample -c FALSE -o SRR3479243
```

## 7. Instruction of Circall simulator

### Introduction
Circall simulator is an R package used to generate both circRNA and tandem RNA sequecing data. 

### Software requirements:
Circall simulator is packed in Circall container.
- R packages version 3.6.0 or latter with following installed packages: GenomicFeatures, Biostrings, and polyester.

### Input parametes:

- BSJ_Info a data frame that contains 5 columns which are: Chr, start_EXONSTART, end_EXONEND, GENEID and cCount. Chr is chromosome name with formated as 1:22, X, Y, Mt. start_EXONSTART is starting position starting exon of circRNA, end_EXONEND is ending position ending exon of circRNA, GENEID is gene ID contains circRNA (used to get gene model) and cCount are number of read pair want to generated for the target circRNA.
- tandem_rate the rate tandem RNA that you wish to simulated eg. 0.05
- error_rate sequencing error rate, defaut is 0.005
- set.seed set seed for reproducibility, defaut is on with 2018
- gtfSqlite path to your annotation file, Sqlite formated (generated by GenomicFeatures)
- genomeFastaFile path to your genome fasta file
- txFastaFile path to your transcript fasta file (cDNA)
- out_name prefix output folders, defaut is "Circall_simuation"
- out_dir the directory contains output, defaut is the current derectory

### Example (assummed that your working directory is the folder contain installed Circall and annotation)
/Users/datn/Documents/DATA2020/PaperNAR/Circall_v0.0.0_linux_x86-64/simulator

Load function into your R
```R
files.sources = list.files("Circall_v0.0.0_linux_x86-64/simulator", full.names = TRUE)
sapply(files.sources, source)

```

Initialize a BSJ_info to genrate CircRNA and tandem RNA
```R
Chr = c(7,7,3,5,17,4,7,1,3,1,17,12,14,10,18,17,5,20,16,17)

start_EXONSTART = c(131113792,99795401,172363413,179296769,36918664,151509200,2188787,51906019,57832924,225239153,76187051,111923075,104490906,101556854,196637,21075331,74981032,60712420,56903641,80730328)

end_EXONEND = c(131128461,99796580,172365904,179315312,36918758,151509336,2270359,51913807,57882659,225528403,76201599,111924628,104493276,101572901,199316,21087123,74998635,60716000,56904648,80772810)

GENEID = c("ENSG00000128585","ENSG00000066923","ENSG00000144959","ENSG00000197226","ENSG00000108294","ENSG00000198589","ENSG00000002822","ENSG00000085832","ENSG00000163681","ENSG00000185842","ENSG00000183077","ENSG00000204842","ENSG00000156414","ENSG00000023839","ENSG00000101557","ENSG00000109016","ENSG00000152359","ENSG00000101182","ENSG00000070915","ENSG00000141556")

cCount = sample(2:20000,20)

BSJ_info = data.frame(Chr = Chr, start_EXONSTART = start_EXONSTART, end_EXONEND = end_EXONEND, GENEID = GENEID, cCount = cCount)
```
Run simulation

```R
simulation = Circall_simulate(BSJ_info = BSJ_info, tandem_rate = 0.05, out_name = "Tutorial", gtfSqlite = "Homo_sapiens.GRCh37.75.sqlite", genomeFastaFile = "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa", txFastaFile = "Homo_sapiens.GRCh37.75.cdna.all.fa", out_dir= "./simulation_test")

```


<!-- ## 6. Circall step by step instruction.

This part is writen for experienced users who would like to to run Circall step by step.
We assume that you have successfully run Circall pipeline in section 5.

### 6.1 Extracting unmmaped reads.

```sh
Circall_v0.0.0_linux_x86-64/linux/bin/Circall_wt -i IndexTranscriptome -1 <(gunzip -c sample_01_1.fasta.gz) -2 <(gunzip -c sample_01_2.fasta.gz) -o outDirWT -p 4
```
The output of this step includes multiple fasta files of unmapped reads corresponding to the number of CPU cores, so, we need to merge them together.

```sh
cat outDirWT/UN_read1_* >  Unmapped_readM_1.fa
cat outDirWT/UN_read2_* >  Unmapped_readM_2.fa

```

### 6.2 Mapping to BSJ database to find BSJ reads.

In this step, we map the extracted unmapped reads againts BSJ databse

```sh
Circall_v0.0.0_linux_x86-64/linux/bin/Circall_bsj -i IndexBSJ -1 Unmapped_readM_1.fa -2 Unmapped_readM_2.fa -o outDirBS -p 4
```

Merge multiple BSJ reads together.

```sh
cat outDirBS/BSJ_read1_* > BSJ_read_merged_1.fa
cat outDirBS/BSJ_read2_* > BSJ_read_merged_2.fa
```

### 6.3 Single-end read filtering.

```sh
Rscript Circall_v0.0.0_linux_x86-64/R/doSingleEndFiltering.R gtfSqlite=Homo_sapiens.GRCh37.75.sqlite CPUNUM=4 outDirBS=outDirBS outFn_SE_filtering_Rdata=Circall_SE_filter_output.RData
```

### 6.4 Generating pseudo transcripts.

```sh
Rscript Circall_v0.0.0_linux_x86-64/R/getPseudoCircRNAsequence.R genomeFastaFile=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa gtfSqlite=Homo_sapiens.GRCh37.75.sqlite txFastaFile=Homo_sapiens.GRCh37.75.cdna.all.fa CPUNUM=4 tandem=TRUE outFn_SE_filtering_Rdata=Circall_SE_filter_output.RData outFn_getPseudoSeq_Rdata=Circall_circRNAinfo_pseudoSeq.RData outFn_getPseudoSeq_fasta=Circall_pseudoSeq.fa
```

### 6.5 Pair-end read filtering.
Index the generated pseudo sequences.
```sh
Circall_v0.0.0_linux_x86-64/linux/bin/TxIndexer -t Circall_pseudoSeq.fa -o pseudo_Idx
```
Map the extracted BSJ reads againts the pseudo sequences index.

```sh
Circall_v0.0.0_linux_x86-64/linux/bin/Circall_pseudo -i pseudo_Idx -1 BSJ_read_merged_1.fa -2 BSJ_read_merged_2.fa -o outDir_pseudo -p 4
```

Merge multiple output files together.

```sh
# note that we need to keep merged_hitheader.txt and merged_mapInfo.txt unchanged, and they must inside the folder of outDir_pseudo to avoid error in the next step. 
cat outDir_pseudo/*.rh > outDir_pseudo/merged_hitheader.txt
rm outDir_pseudo/*.rh
cat outDir_pseudo/mapInfo_*.txt > outDir_pseudo/merged_mapInfo.txt
rm outDir_pseudo/mapInfo_*.txt
```

Run Rscript to performe pair-end read filtering

```sh
Rscript Circall_v0.0.0_linux_x86-64/R/doPairEndFiltering.R outDirWT=outDirWT outDir_pseudo=outDir_pseudo outFn_SE_filtering_Rdata=Circall_SE_filter_output.RData outFn_PE_filtering_Rdata=Circall_PE_filter_output.RData
```
Now, your circRNA candidate list is in Circall_PE_filter_output.RData. If you would like to stop without fdr2d (your input data is Rnase R treated sample), you need to export results to a txt file.
```sh
Rscript Circall_v0.0.0_linux_x86-64/R/export.R outFn_PE_filtering_Rdata=Circall_PE_filter_output.RData outFn_circRNA_final=Circall_final.txt
```

### 6.6 Run fdr2d.

```sh
Rscript Circall_v0.0.0_linux_x86-64/R/getFdr.R outFn_PE_filtering_Rdata=Circall_PE_filter_output.RData depDataFile=Circall_v0.0.0_linux_x86-64/Data/Circall_depdata_human.RData outFn_circRNA_final=Circall_final.txt
```


 -->
