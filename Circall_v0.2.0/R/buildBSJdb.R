#!/usr/bin/env Rscript

# USAGE: this script is used to buil peseudo BSJ sequences.
# Input:
# genomeFastaFile: genome sequence
# gtfSqlite: sqlite of the gtf file of gene annotation
# chromRef: list of selected chromosomes separated by comma ","
# bsjDist: Distance to the the back splicing junction of a circular RNA (default setting: bsjDist=250)
# maxReadLen: maximum read length of RNA sequencing (default setting: maxReadLen=150)
# output: output file name
# 
# Output

# systax: Rscript buildBSJdb.R gtfSqlite=$gtfSqlite genomeFastaFile=$genomeFastaFile bsjDist=$bsjDist maxReadLen=$maxReadLen output=$output

###
genomeFastaFile=gtfSqlite=maxReadLen=bsjDist=out_fasta=chromRef=NA
args = commandArgs(trailingOnly=TRUE)
#cat("\nNumber of arguments: ",length(args))
#cat("\nList of arguments: ",args)

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="genomeFastaFile") genomeFastaFile=res[2]
	if (res[1]=="gtfSqlite") gtfSqlite=res[2]
  if (res[1]=="chromRef") chromRef=res[2]
  if (res[1]=="maxReadLen") maxReadLen=res[2]
	if (res[1]=="bsjDist") bsjDist=res[2]
	if (res[1]=="output") out_fasta=res[2]
}



#check input information
validatedCommand=TRUE
if (is.na(genomeFastaFile)){
  cat("\nThere is no genome sequence. Stop!")
  validatedCommand=FALSE
}
if (is.na(gtfSqlite)){
  cat("\nThere is no sqlite file. Stop!")
  validatedCommand=FALSE
}
if (is.na(out_fasta)){
  cat("\n There is no output file name, use the default name.")
  prefix="Circall"
  #out_exon_intron = paste0(prefix,"_intronexonSeq.RData")
  out_fasta = paste0(prefix,"_BSJ_sequences.fa")
  #out_Rdata = paste0(prefix,"_BSJ_names.Rdata")  
}
if (is.na(chromRef)){
  cat("\n There is no setting for the list of selected chromosomes, use all chromosomes")
}else{
  #chromRef="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"
#  cat("\n Selected chromosomes: ",chromRef)
  chromRef=trimws(unlist(strsplit(chromRef,",")))
}

if (is.na(bsjDist)){
  bsjDist=250
  cat("\n There is no setting for the distance to the the back splicing junction, use the default setting (bsjDist=250).")
}
if (is.na(maxReadLen)){
  maxReadLen=150
  cat("\n There is no setting for the maximum read length of RNA sequencing,  use the default setting (maxReadLen=150).")
}

if (!validatedCommand) stop("\nStop without results!!!")


if (validatedCommand){
  cat("\n-----")
  cat("\nParameter settings:")
  cat("\n genomeFastaFile=",genomeFastaFile)
  cat("\n gtfSqlite=",gtfSqlite)
  cat("\n chromRef=",paste(chromRef,collapse=","))   
  cat("\n bsjDist=",bsjDist)
  cat("\n maxReadLen=",maxReadLen)
  cat("\n output=",out_fasta)
}

# some variable must be integer
maxReadLen=as.integer(maxReadLen)
bsjDist= as.integer(bsjDist)

# output names
#prefix="Circall"
#out_exon_intron = paste0(prefix,"_intronexonSeq.RData")
#out_fasta = paste0(prefix,"_BSJ_sequences.fa")
#out_Rdata = paste0(prefix,"_BSJ_names.Rdata")

# #load packages
# suppressMessages(suppressWarnings(require(foreach)))
# suppressMessages(suppressWarnings(require(doParallel)))
# registerDoParallel(cores=CPUSNUM) #using param setting CPUs number
suppressMessages(suppressWarnings(require(Biostrings)))
suppressMessages(suppressWarnings(require(GenomicFeatures)))


anntxdb <- loadDb(gtfSqlite)

fasta = readDNAStringSet(genomeFastaFile)
fasta_chrnames=sapply(names(fasta), function(x) unlist(strsplit(x," "))[1])

intron.info=list()
exon.info=list()
for (chrID in 1:length(fasta_chrnames)){
  myex=exons(anntxdb,columns=c("EXONID"),filter=list(exon_chrom=fasta_chrnames[chrID]))
  if (length(myex)>0){
    startPos=unique(start(myex))
    endPos=unique(end(myex))
    plusStrand=strand(myex)=="+"
    
    intronS=substring(fasta[chrID],startPos-bsjDist,startPos-1)
    intronE=substring(fasta[chrID],endPos+1,endPos+bsjDist)
    intronSeq.start=intronS
    intronSeq.end=intronE
    exonStrand=plusStrand
    intron.info[[fasta_chrnames[chrID]]]=list(intronSeq.start=intronSeq.start,intronSeq.end=intronSeq.end,exonID=elementMetadata(myex)$EXONID,exonStrand=exonStrand, startPos=startPos,endPos=endPos)
    
    exonS=substring(fasta[chrID],startPos,startPos+bsjDist-1)
    exonE=substring(fasta[chrID],endPos-bsjDist+1,endPos)
    exonSeq.start=exonS
    exonSeq.end=exonE
    exonStrand=plusStrand
    exon.info[[fasta_chrnames[chrID]]]=list(exonSeq.start=exonSeq.start,exonSeq.end=exonSeq.end,exonID=elementMetadata(myex)$EXONID,exonStrand=exonStrand,startPos=startPos,endPos=endPos)
#    cat("DONE processing chromosome ",fasta_chrnames[chrID],"\n")
  }
}

############################################################
# Finished getting exonintronSeq
# Begin detecting canonical splicing
############################################################

convertReverseComplement<-function(DNAseq){
  DNAarr=unlist(strsplit(DNAseq,""))
  #reverse
  DNAarr=rev(DNAarr)
  #complement
  Aid=which(DNAarr=="A")
  Tid=which(DNAarr=="T")
  Gid=which(DNAarr=="G")
  Cid=which(DNAarr=="C")
  DNAarr[Aid]="T"
  DNAarr[Tid]="A"
  DNAarr[Gid]="C"
  DNAarr[Cid]="G"
  #result
  DNAseqRc=paste(DNAarr,collapse = "")
  return(DNAseqRc) 
}


genes.all=genes(anntxdb)
genes.exon.all=select(anntxdb, keys=genes.all$gene_id, columns=c("GENEID","EXONID","EXONSTART","EXONEND","EXONSTRAND","EXONCHROM"), keytype = "GENEID")

if (is.na(chromRef)){ #use all chromosomes
  chromRef=names(exon.info)
}
#chrVec=names(exon.info)
#chrVec=chrVec[1:25] #chr1-22,X,Y and MT
chrVec=chromRef

bs.names.vec=NULL
for (chr in chrVec){
  cat("\n Begin processing Chromosome: ",chr)
  exon.chr=exon.info[[chr]]
  genes.chr=genes(anntxdb, filter=list(tx_chrom = chr))
  cat("\n # of genes ",length(genes.chr),"\n")
  genes.exon.map=select(anntxdb, keys=names(genes.chr), columns=c("GENEID","EXONID","EXONSTART","EXONEND","EXONSTRAND","EXONCHROM"), keytype = "GENEID")
  #remove exons from different chromosomes: some genee have exons from two chromosomes, e.g. AKAP17A
  genes.exon.map=genes.exon.map[genes.exon.map$EXONCHROM==chr,]
  cat("# of exons ",nrow(genes.exon.map),"\n")
  #get intron information
  intron.chr=intron.info[[chr]]
  #get two bases from introns of the startPos and endPos
  prime3seq=sapply(intron.chr$intronSeq.start, function(x) substring(x,nchar(x)-1))
  prime5seq=sapply(intron.chr$intronSeq.end, function(x) substring(x,1,2))
  names(prime5seq)=NULL
  names(prime3seq)=NULL  
  #get reversed complement sequence of the intron: 5 prime and 3 prime are exchanged
  prime5seq.rc=sapply(intron.chr$intronSeq.start, function(x) substring(convertReverseComplement(x),1,2))
  prime3seq.rc=sapply(intron.chr$intronSeq.end, function(x) substring(convertReverseComplement(x),nchar(x)-1))
  names(prime5seq.rc)=NULL
  names(prime3seq.rc)=NULL
  for (g in names(genes.chr)){
    #g=names(genes.chr)[1]
    GEmat=genes.exon.map[genes.exon.map$GENEID==g,]
    GEmat=GEmat[order(GEmat$EXONID),] #sort exons by increasing order for forward strand
    #get exon info
    pick=exon.chr$exonID %in% GEmat$EXONID

    GEmat$EXONLEN=GEmat$EXONEND-GEmat$EXONSTART+1
    s=unique(GEmat$EXONSTART)
    e=unique(GEmat$EXONEND)
    #select the exonID of the smallest exon length
    exID.s=sapply(s,function(x){pick=GEmat$EXONSTART==x;return(GEmat$EXONID[pick][which.min(GEmat$EXONLEN[pick])])})
    exID.e=sapply(e,function(x){pick=GEmat$EXONEND==x; return(GEmat$EXONID[pick][which.min(GEmat$EXONLEN[pick])]) }  )
    #
    exIDMat=expand.grid(exID.s,exID.e)      
    itvMat=expand.grid(s,e)
    #number of exons in the regions of s-e, this similar as itvMat$Var2-itvMat$Var1 >=0. The combination is valid if EXONEND > EXONSTART
    itvExnum=apply(itvMat,1, function(x) sum(x[1] <= GEmat$EXONSTART & x[2] >= GEmat$EXONEND))
    valid.itv=which(itvExnum > 0)
    ### get matched posistions
    # sequences from exon.chr$exonSeq.start corresponding to exon.chr$startPos, similarly to endPos
    pick.s=match(itvMat$Var1[valid.itv],exon.chr$startPos)
    pick.e=match(itvMat$Var2[valid.itv],exon.chr$endPos)
    #CHECK canonical criteria: if endPos(GT) and startPos(AG)
    #Statistics of splice pairs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC113136/
    if (GEmat$EXONSTRAND[1]=="+"){
      #isCS=prime3seq[pick.s]=="AG" & prime5seq[pick.e]=="GT" #GT-AG with 99.24%
      isCS.GT_AG=prime3seq[pick.s]=="AG" & prime5seq[pick.e]=="GT" #GT-AG with 99.24%
      isCS.GC_AG=prime3seq[pick.s]=="AG" & prime5seq[pick.e]=="GC" #GC-AG with 0.69%
      isCS.AT_AC=prime3seq[pick.s]=="AC" & prime5seq[pick.e]=="AT" #AT-AC with 0.05%
      isCS= isCS.GT_AG | isCS.GC_AG | isCS.AT_AC
    }
    #if minus strand - use reverse complement sequence
    if (GEmat$EXONSTRAND[1]=="-"){ 
      #isCS=prime3seq.rc[pick.e]=="AG" & prime5seq.rc[pick.s]=="GT"
      isCS.GT_AG=prime3seq.rc[pick.e]=="AG" & prime5seq.rc[pick.s]=="GT" #GT-AG with 99.24%
      isCS.GC_AG=prime3seq.rc[pick.e]=="AG" & prime5seq.rc[pick.s]=="GC" #GC-AG with 0.69%
      isCS.AT_AC=prime3seq.rc[pick.e]=="AC" & prime5seq.rc[pick.s]=="AT" #AT-AC with 0.05%
      isCS= isCS.GT_AG | isCS.GC_AG | isCS.AT_AC
    }
    if (sum(isCS)>0) {
      ###create back splicing sequence: seq-end of end-exon + seq-start of start
      #if (GEmat$EXONSTRAND[1]=="+") 
      bs.seq=paste(substr(exon.chr$exonSeq.end[pick.e],bsjDist-maxReadLen+2,bsjDist),substr(exon.chr$exonSeq.start[pick.s],1,maxReadLen-1),sep="")
      #if minus strand: still keep the same same as plus strand, but the read must be mapped in RC
      #if (GEmat$EXONSTRAND[1]=="-") bs.seq=paste(substr(exon.chr$exonSeq.end[pick.e],bsjDist-maxReadLen+2,bsjDist),substr(exon.chr$exonSeq.start[pick.s],1,maxReadLen-1),sep="")
      #Txname: chr_startPos_endPos_gene_starExonID_endExonID tx
      bs.names=paste(">",chr,"__",itvMat$Var1[valid.itv],"__",itvMat$Var2[valid.itv],"__",g,"__",exIDMat$Var1[valid.itv],"__",exIDMat$Var2[valid.itv]," bs",sep="")
      #keep only the canonical splicing
      bs.seq=bs.seq[isCS]
      bs.names=bs.names[isCS]
      bs.final=c(bs.names,bs.seq)
      bs.final[seq(1,length(bs.final),2)]=bs.names
      bs.final[-seq(1,length(bs.final),2)]=bs.seq      
      write.table(bs.final, file=out_fasta, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
      bs.names.vec=c(bs.names.vec,bs.names)
    }
  }
}

#cat("Exporting exonintronSeq Rdata file! \n")
#save(exon.info,intron.info, file=out_exon_intron)
#
#cat("\n Exporting BSJ_name Rdata file! \n")
#save(bs.names.vec, file=out_Rdata)

rm(list=ls())

