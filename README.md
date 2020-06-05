
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
- dep -- data contain depleted circRNAs: specify location of deletep data that is used as training data for estimation of fdr2d.
- p -- number of thread: Defaut is 4
- tag -- tag name of results:  Defaut is "Sample".
- td -- generation of tandem sequences: TRUE/FALSE value, defaut is TRUE
- c -- clean intermediate data: TRUE/FALSE value, defaut is TRUE
- o -- output folder: Defaut is current directory

## 6. Circall step by step instruction.

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

### 6. Pair-end read filtering.
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
Rscript Circall_v0.0.0_linux_x86-64/R/doPairEndFiltering.R outDirWT=outDirWT outDir_pseudo=outDir_pseudo outFn_SE_filtering_Rdata=Circall_SE_filter_output.RData outFn_PE_filtering_Rdata=Circall_PE_filter_output.RData CPUNUM=4
```
Now, your circRNA candidate list is in Circall_PE_filter_output.RData

