#load arguments
CPUNUM=NA
args = commandArgs(trailingOnly=TRUE)

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="genomeFastaFile") genomeFastaFile=res[2]
	if (res[1]=="gtfSqlite") gtfSqlite=res[2]
	if (res[1]=="txFastaFile") txFastaFile=res[2]
	if (res[1]=="CPUNUM") CPUNUM=res[2]

	if (res[1]=="outDirBS") outDirBS=res[2]
	if (res[1]=="outFn_SE_filtering_Rdata") outFn_SE_filtering_Rdata=res[2]
	
	if (res[1]=="tandem") tandem=res[2]
	if (res[1]=="outFn_getPseudoSeq_Rdata") outFn_getPseudoSeq_Rdata=res[2]
	if (res[1]=="outFn_getPseudoSeq_fasta") outFn_getPseudoSeq_fasta=res[2]

	if (res[1]=="outDirWT") outDirWT=res[2]
	if (res[1]=="outDir_pseudo") outDir_pseudo=res[2]
	if (res[1]=="outFn_PE_filtering_Rdata") outFn_PE_filtering_Rdata=res[2]

	if (res[1]=="depDataFile") depDataFile=res[2]
	if (res[1]=="outFn_circRNA_final") outFn_circRNA_final=res[2]		

}


if (is.na(CPUNUM)) CPUNUM=1
CPUNUM = as.integer(CPUNUM)

#load packages
suppressMessages(suppressWarnings(require(foreach)))
suppressMessages(suppressWarnings(require(doParallel)))
suppressMessages(suppressWarnings(require(Biostrings)))
suppressMessages(suppressWarnings(require(GenomicFeatures)))
#cores
registerDoParallel(cores=CPUNUM)
#load R functions
source("/path/to/R/Circall_functions.R")