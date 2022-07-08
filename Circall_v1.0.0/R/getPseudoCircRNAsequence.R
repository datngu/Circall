# This module used to get pseudo circRNA sequences from results of SE filtering

# Input: outFn_SE_filtering_Rdata, gtfSqlite annotation, genome fasta file and transcript (cDNA) fasta file.
# Output: fasta and Rdata that contains pseudo sequences.

# Syntax doSingleEndFiltering.R outDirBS=path/to/outDirBS outFn_SE_filtering=path/to/SE_filtering.Rdata CPUNUM=cpu_number genomeFastaFile=path/to/genome.fasta gtfSqlite=path/to/gtf.sqlite txFastaFile=path/to/transcript_cDNA.fa outFn_getPseudoSeq_Rdata=PseudoSeq.Rdata outFn_getPseudoSeq_fasta=PseudoSeq.fa


# UPDATE NOTES:
# 03 Dec 2019: new format of tandem name, sep stat, break, end by "--" in order to avoid similar with tx name of rice
# 19 Apr 2021/Dat: fix bug on hg38 annotation, tx names with '.version' after tx name

#########################
rm(list = ls())
# process arguments and load functions
source("/path/to/R/load_packages.R", print.eval = FALSE)

# load annotation
# genes.exon.all=process_Sqlite(gtfSqlite)
 anntxdb <- loadDb(gtfSqlite)
 genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
 genes.exon.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID", "TXNAME", "EXONID", "EXONSTART", "EXONEND", "EXONSTRAND", "EXONCHROM"), keytype = "GENEID")))
 genes.exon.all$LENGTH = genes.exon.all$EXONEND - genes.exon.all$EXONSTART + 1
# load genome fasta
# fasta_genome=process_genome(genomeFastaFile)
 fasta_genome = readDNAStringSet(genomeFastaFile)
 chnames = sapply(names(fasta_genome), function(x) unlist(strsplit(x, " "))[1])
 names(fasta_genome) = chnames
 
#load tx.all.fasta and extract txname to dataframe
tx.all.fasta = readDNAStringSet(txFastaFile)
# to get seqs from cDNA
# 19Apr2021/Dat: fix bug on hg38 annotation, tx names with '.version' after tx name
tx.all.NAME = sapply(names(tx.all.fasta),function(x) unlist(strsplit(x," "))[1])
tx.all.NAME = sapply(tx.all.NAME,function(x) unlist(strsplit(x,".", fixed = T))[1])

SE_filtering_res = load_SE_filtering(outFn_SE_filtering_Rdata)
readLen = SE_filtering_res$fragInfo$readlen

# collecting circRNA list from SE_filtering
circRNAs_candidates = as.character(SE_filtering_res$circRNA.names)
circInfo = get_circInfo(circRNAs_candidates)
# get pseudo sequences

#par.pseudoSeq.fa = foreach(i = 1:nrow(circInfo),.combine = rbind) %dopar% {
#  #par.pseudoSeq.fa=foreach(i = 1:20,.combine=rbind) %dopar% {
#  mycirR = circInfo[i,]
#  mygene = as.character(mycirR$GENEID)
#  myexonInfo = genes.exon.all[genes.exon.all$GENEID == mygene, ]
#  start = mycirR$start_EXONSTART
#  end = mycirR$end_EXONEND
#  circRNA_name = circRNAs_candidates[i]
#
#  if(tandem == "TRUE"){
#    res1 = getCircPseudoSequence_cDNA(myexonInfo, circRNA_name, start, end, num = -1, readLen = readLen)
#    res2 = getTandemSequence_cDNA(myexonInfo, circRNA_name, start, end, num = -1, readLen = readLen)
#    res = rbind(res1, res2)
#  }else{
#    res = getCircPseudoSequence_cDNA(myexonInfo, circRNA_name, start, end, num = -1, readLen = readLen)
#  }
#  return(res)
#}
#
#
#pseudoSeq.fa = DNAStringSet("ATCG") #start with a dumpy DNA sequence
#for (i in 1:nrow(par.pseudoSeq.fa)){
#  circRNA.fa = DNAStringSet(par.pseudoSeq.fa[i,1])
#  names(circRNA.fa) = par.pseudoSeq.fa[i,2]
#  pseudoSeq.fa =c (pseudoSeq.fa, circRNA.fa)
#}
#pseudoSeq.fa = pseudoSeq.fa[-1] #remove the dumpy DNA sequence

par.pseudoSeq.fa=foreach(i = 1:nrow(circInfo)) %dopar% {
# par.pseudoSeq.fa=lapply(1:nrow(circInfo), function(i){
  #par.pseudoSeq.fa=foreach(i = 1:20,.combine=rbind) %dopar% {
  mycirR=circInfo[i,]
  mygene=as.character(mycirR$GENEID)
  myexonInfo=genes.exon.all[genes.exon.all$GENEID==mygene,]
  start=mycirR$start_EXONSTART
  end=mycirR$end_EXONEND
  circRNA_name=circRNAs_candidates[i]

  if(tandem == "TRUE"){
    res1=getCircPseudoSequence_cDNA(myexonInfo, circRNA_name, start, end, num=-1,readLen=readLen)
    res2=getTandemSequence_cDNA(myexonInfo, circRNA_name, start, end, num=-1,readLen=readLen)
    #res3=getTandemSequence(myexonInfo, circRNA_name, start, end, num=-1)
    fa = c(res1$fa, res2$fa)
    name = c(res1$name, res2$name)
    res= list("fa" = fa, "name" = name)
  }else{
    res=getCircPseudoSequence_cDNA(myexonInfo, circRNA_name, start, end, num=-1,readLen=readLen)
  }
  return(res)
}

pseudoSeq.fa = lapply(par.pseudoSeq.fa, function(p) p$fa)
pseudoSeq.name = lapply(par.pseudoSeq.fa, function(p) p$name)

pseudoSeq.fa = unlist(pseudoSeq.fa)
pseudoSeq.name = unlist(pseudoSeq.name)

pseudoSeq.fa = DNAStringSet(pseudoSeq.fa)
names(pseudoSeq.fa) = pseudoSeq.name

# export to file
save(pseudoSeq.fa, circInfo, file = outFn_getPseudoSeq_Rdata)
# export fasta of pseudoSeq to fasta file
writeXStringSet(pseudoSeq.fa, outFn_getPseudoSeq_fasta)

# #export fasta of the ONLY circRNAs (CIRC) to file
# fastaFile_onlyCIRC=paste(tagStr,"_onlycircRNA.fa",sep="")
# pick=grepl("--",names(pseudoSeq.fa))
# onlyCIRCFasta=pseudoSeq.fa[!pick]
# writeXStringSet(onlyCIRCFasta, fastaFile_onlyCIRC)

rm(list=ls())