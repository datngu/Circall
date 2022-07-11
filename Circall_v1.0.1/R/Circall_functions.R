# 22 Nov 2018: fix the prolem of changing matrix format when select a single row/column
# 25 Otc 2018: new format of normID: chr__start__end
# 25 Otc 2019: update read name for PE filtering "___??" - Dat
# 03 Dec 2019: new format of tandem name, sep stat, break, end by "--" in order to avoid similar with tx name of rice - Dat
# 06 Dec 2019: update filtering phase 1, for 1 exon circRNA: min hit lenght is x2 the circRNA length before finding min("readlen","f_reflen","circ.min.len"). - Dat
# add median.len, number of tx, etc.. to out put, and export output as txt file.- Dat
# 04 Jun 2020: speed up the process of extracting sequences for circRNA and tandem - Thinh
# 19 Apr 2021/Dat: fix bug on hg38 annotation, tx names with '.version' after tx name
# 05 July 2022: add data.table for faster reading.
# contact: Trung Nghia Vu (TrungNghiaVu@ki.se)

SE_filtering <- function(inPath, gtfSqlite, outFn = NULL, read_mapProp = 0.99, read_anchorLen = 10, eqc_txMinProp = 1, eqc_minCount = 1, hard_hitlen = 50){

  fileList = list.files(path = inPath, pattern = "UN_read", full.names = TRUE)
  fragInfo = read.csv(paste(inPath, "/fragmentInfo.txt", sep = ""), sep = "\t")
  # input read data
  read_dat=NULL
  for (fn in fileList){
  a_read_dat = fread(fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  a_read_dat = as.data.frame(a_read_dat)
  read_dat = rbind(read_dat, a_read_dat)
  }
  colnames(read_dat) = c("Header", "Direction", "QueryPos", "HitLen", "HitPos", "Txname", "Reflen")
  read_dat$dis5 = read_dat$Reflen/2 - read_dat$HitPos
  read_dat$dis3 = read_dat$HitPos + read_dat$HitLen - read_dat$Reflen/2

  #input data to get circLen from read_dat
  all_circInfo = get_circInfo(unique(as.character(read_dat$Txname)))
  genes.exon.all = process_Sqlite(gtfSqlite)
  all_circLen = get_circLen(genes.exon.all, all_circInfo)
  #processing of read_dat
  read_dat$normID = sapply(read_dat$Txname, function(x){
    y = gregexpr("__", x)
    z = y[[1]][3]-1
    return(substring(x, 1, z))
  })
  
  match_idx = match(read_dat$normID, all_circLen$normID)
  my_dat = all_circLen[match_idx, ]
  read_dat = cbind(read_dat, my_dat)
  read_dat$min.circlen = ifelse(read_dat$max_exon_num == 1, read_dat$min.circlen*2, read_dat$min.circlen)
  read_dat$readlen = fragInfo$readlen
  read_dat$f_reflen = read_dat$Reflen/2 + read_anchorLen
  read_dat$min_hitlength = read_mapProp*apply(read_dat[ ,c("readlen", "f_reflen", "min.circlen")], 1, min)
  read_dat$min_hitlength = ifelse(read_dat$min_hitlength >= hard_hitlen, read_dat$min_hitlength, hard_hitlen)
  #add hard hitLen to limit for too short circRNAs 

  #get equivalence classes
  eqc_dat = read.table(paste(inPath, "/rawCount.txt", sep=""), sep="\t", header=TRUE)
  #save raw data
  read_raw = read_dat
  eqc_raw = eqc_dat

  #filter read_dat
  read_dat = read_dat[read_dat$HitLen >= read_dat$min_hitlength, ] # reads that mapped read_mapProp of read-length
  pick = read_dat$dis3 > read_anchorLen & read_dat$dis5 > read_anchorLen # reads must overlap more than read_anchorLen bases of the left (right) sequence at the break-point  
  read_dat = read_dat[pick, ]

  #do some simple filters for equivalence-class data
  pick = eqc_dat$Weight >= eqc_txMinProp & eqc_dat$Count >= eqc_minCount # reads of the eqc must map uniquely to the transcript; and eqc must contain least 1 reads 
  eqc_dat = eqc_dat[pick, ]
  # get overlapping with the results from the read_dat
  pick = eqc_dat$Transcript %in% unique(read_dat$Txname)
  eqc_dat = eqc_dat[pick, ]
  #dim(eqc_dat)
  #rank by expression
  eqc_dat = eqc_dat[order(eqc_dat$Count, decreasing = TRUE), ]
  eqc_dat$normID = sapply(eqc_dat$Transcript, function(x){
    y = gregexpr("__", x)
    z = y[[1]][3]-1
    return(substring(x, 1, z))
  })
  
  eqc_dat$median.circlen = all_circLen$median.circlen[match(eqc_dat$normID, all_circLen$normID)]


  match_idx = match(eqc_dat$normID, all_circLen$normID)
  my_dat = all_circLen[match_idx, ]
  eqc_dat = cbind(eqc_dat, my_dat)

  fragment_count = table(read_dat$Txname)
  eqc_dat$fragment_count_SE = fragment_count[match(eqc_dat$Transcript, names(fragment_count))]
  eqc_dat$fragment_count_SE = as.integer(eqc_dat$fragment_count)
  #writing file txt for quick look
  #write.table(eqc_dat, file = "Phase1Filterd_circRNA.txt", sep = "\t",row.names = FALSE)
  circRNA.names = unique(as.character(eqc_dat$Transcript))
  #length(circRNA.names)
  #save to file
  if(!is.null(outFn)) save(circRNA.names, eqc_dat, read_dat, fragInfo, read_raw, eqc_raw, file = outFn)
  res = list(circRNA.names = circRNA.names, eqc_dat = eqc_dat, read_dat = read_dat, fragInfo = fragInfo, read_raw = read_raw, eqc_raw = eqc_raw)
  
  return(res)
}


#load phase 1 - information
loadPhase1Info <- function(dataname){
  load(dataname)
  res = list(circInfo = circInfo,circRNA.tx.fa = circRNA.tx.fa)
  return(res)
}

load_SE_filtering <- function(dataname){
  load(dataname)
  #res=list(circRNA.names=circRNA.names,eqc_dat=eqc_dat,read_dat=read_dat,fragInfo=fragInfo)
  res = list(circRNA.names = circRNA.names, eqc_dat = eqc_dat, read_dat = read_dat, fragInfo = fragInfo, read_raw = read_raw, eqc_raw = eqc_raw)
  return(res)
}





# PE filter function
#Dat 15/08/2019


load_read2_2_file <- function(path_to_header, path_to_mapInfo){
  Txname = fread(path_to_mapInfo, sep = "\t", stringsAsFactors = FALSE, header=FALSE)
  Header = fread(path_to_header, sep = "\t", stringsAsFactors = FALSE, header=FALSE)
  read2 = cbind(Header, Txname)
  read2 = as.data.frame(read2)
  colnames(read2) = c("Header", "Txname")
  return(read2)
}

load_resPhase1 <- function(path_to_Res){
  load(path_to_Res)
  res = list(read_raw = read_raw, eqc_raw = eqc_raw, eqc_dat = eqc_dat)
  return(res)
}


#10-Feb-2020 update tandem filtering
PE_filter <- function(PEread, resPhase1, fragInfo_wt, read_mapProp = 0, read_anchorLen = 10, eqc_txMinProp = 1, min_fragment_count = 1){
# #input for testing
  # PEread=read2.minus
  # resPhase1=res.minus
  # read_mapProp=0
  # read_anchorLen=10
  # eqc_txMinProp=1
  # min_fragment_count=2
  # readlen=100
# ########
  #PEread=PEread2


  #process single-end mapping (SE)
  SE = resPhase1$read_raw
  SE$Header = gsub("___11|___12|___21|___22", "", SE$Header)
  SE$Txname = as.character(SE$Txname)
  #SE=SE[SE$HitLen>= SE$min_hitlength,] #apply stringent like phase 1 or not?
  pick = SE$dis3 > read_anchorLen & SE$dis5 > read_anchorLen
  SE = SE[pick, ]
  SE$header_TxnameSE = paste(SE$Header, SE$Txname,sep = "-") #combine read headers and mapped transcript name for comparison with paired-end mapping

  #process paired-end (PE) mapping information
  PEread$Header = gsub("___11|___12|___21|___22", "", PEread$Header)
  #this morment, one Txname may contains more than one Txname due to multiple hits. e.g one gene model produces many isoforms both tandem and circRNAs => these following step to isolate them out. 
  
  PEread = PEread[PEread$Header %in% SE$Header, ] #pre filter out unneeded reads to speed up!

  # separate reads mapped unique to tandem
  pick = grepl("++", PEread$Txname, fixed = TRUE)
  td_read = PEread[!pick, ]
  ################

  #strsplit out mapped transcripts
  PE.tx = sapply(PEread$Txname, strsplit, split=' ')
  names(PE.tx) = PEread$Header
  TxnamePE = unlist(PE.tx, use.names = FALSE)
  names(TxnamePE) = rep(names(PE.tx), lengths(PE.tx)) #length() and lengths() are different
  #now, one row is stand for 1 mapped transcripts. Read headers repeated multiple times!
  PE = as.data.frame(TxnamePE)
  PE$Header = names(TxnamePE)
  #keep orignal transcript PE names
  PE$TxnamePE0 = PE$TxnamePE
  #strim out tails of transcript PE names and then combine with read headers to compare with single end mapping
  PE$TxnamePE = sapply(PE$TxnamePE, function(x){
    y = gregexpr("__", x)
    #s=y[[1]][3]+2
    e = y[[1]][6] - 1
    return(substring(x, 1, e))
  }) 

  PE$header_TxnamePE = paste(PE$Header, PE$TxnamePE, sep = "-")
  #re-order columns
  PE = PE[c("Header", "header_TxnamePE", "TxnamePE0", "TxnamePE")]
  PE_raw = PE # back up PE_raw
  pick = PE$header_TxnamePE %in% SE$header_TxnameSE #merging information bw SE and PE mapping
  PE = PE[pick,]
  PE$HitLen_SE = SE$HitLen[match(PE$header_TxnamePE, SE$header_TxnameSE)]
  

  #PE.BK = PE
  #PE = PE.BK
  #remove duplicated reads
  PE = PE[!duplicated(PE$header_TxnamePE), ] 

  #tandem_read=get_unique_tandem_read(PE)
  pick = PE$Header %in% td_read$Header
  #pick out tandem unique mapped read
  PE_td = PE[pick,]
  PE = PE[!pick,]
  #now doing recount for eqc_dat
  eqc_dat = resPhase1$eqc_dat # using eqc_matched with median.len to take circRNA length infomation

  tandem_unique_count = table(PE_td$TxnamePE)
  fragment_count = table(PE$TxnamePE)

  eqc_dat$fragment_count = fragment_count[match(eqc_dat$Transcript, names(fragment_count))]
  eqc_dat$fragment_count = as.integer(eqc_dat$fragment_count)

  eqc_dat$tandem_unique_count = tandem_unique_count[match(eqc_dat$Transcript, names(tandem_unique_count))]
  eqc_dat$tandem_unique_count = as.integer(eqc_dat$tandem_unique_count)
  eqc_dat$tandem_unique_count[is.na(eqc_dat$tandem_unique_count)] = 0

  ### get overlapping with the results from the read_dat
  pick = eqc_dat$Transcript %in% unique(PE$TxnamePE)
  eqc_dat = eqc_dat[pick, ]
  ### do some simple filters for equivalence-class data
  pick = eqc_dat$Weight >= eqc_txMinProp & eqc_dat$fragment_count >= min_fragment_count
  eqc_dat = eqc_dat[pick, ]
  #exclude candidate with tandem_unique_count>= 1/2 fragment_count
  pick = eqc_dat$fragment_count >= 3*eqc_dat$tandem_unique_count
  eqc_dat = eqc_dat[pick, ]
  #make normID to easy compare results
 #  eqc_dat$normID =sapply(eqc_dat$Transcript,function(x){
 #    y=gregexpr("__",x)
 #    z=y[[1]][3]-1
 #    return(substring(x,1,z))
 #  })
  circRNAs_candidates = eqc_dat$Transcript

  eqc_dat$fragment_count_norm=eqc_dat$fragment_count/fragInfo_wt$numObservedFragments*10^6

  res = list(circRNAs_candidates = circRNAs_candidates, fragInfo_wt=fragInfo_wt, eqc_dat = eqc_dat, PE = PE, PE_raw = PE_raw)
  return(res)
}


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


getCircPseudoSequence_cDNA <-function (myexonInfo, circRNA_name, start, end, num=-1,readLen){
  #output
  circRNA.tx.fa = c() #start with a dumpy DNA sequence - serial mode
  circRNA.tx.name = c() #start with a dumpy DNA sequence - serial mode
  #get all exons completely inside the boundary
  pick1=myexonInfo$EXONSTART == start
  pick2=myexonInfo$EXONEND == end
  mytx=intersect(unique(myexonInfo[pick1,]$TXNAME), unique(myexonInfo[pick2,]$TXNAME))
  #find all transcripts and generate circRNA from those transcripts
  txNum=length(mytx)

  if (txNum > 0){ # if existing transcripts containing both EXONSTART and EXONEND
    #atx=sample(mytx,1)
    if (num > 0){ #if fixed the number of selected tx, normally num=1 for the simulator tool
      num=min(num,length(mytx))
      mytx=sample(mytx,num)
    }
    
    for (atx in mytx){ #do for every tx
      mytxex=myexonInfo[myexonInfo$TXNAME==atx,]
      mytxex=mytxex[order(mytxex$EXONSTART),]
      mytxex.len = mytxex$LENGTH

      #covert to transcript location
      mytxex.end.all =cumsum(mytxex$LENGTH)
      mytxex.stat.all = mytxex.end.all - mytxex.len +1
      #select the coding area inside the circRNA region
      pick= which(mytxex$EXONSTART >= start & mytxex$EXONEND <= end) #pick by indexes of exons

      startR= mytxex.stat.all[pick]
      endR= mytxex.end.all[pick]
      ## generate pseudo circularRNAs sequences
      # - strand
      if( mytxex$EXONSTRAND[1] == "-"){
        #pick transcript sequence
        transcript = convertReverseComplement(as.character(tx.all.fasta[which(tx.all.NAME == atx)]))
      }else{
        transcript = as.character(tx.all.fasta[which(tx.all.NAME == atx)])
      }
      
      x=substring(transcript,startR,endR)
      txSeq=paste(x,collapse="")

      # new create a pseudo-circularRNA sequence by adding l-1 nucleotides from the end of the last exon to the begin of the first exon
      l_1=readLen -1 -1
      if (l_1 > nchar(txSeq)) l_1=nchar(txSeq) - 1 #if the circRNA length < l_1
      circRNA.seq=paste(substring(txSeq, nchar(txSeq)-l_1),txSeq,sep="")
      circRNA.as_name=paste(atx,"++",sep="") #do not forget to put the "circRNA", otherwises the format of the output file rawcount.txt of FuSeq is broken  
      #add to result
      circRNA.tx.fa=c(circRNA.tx.fa,circRNA.seq)
      circRNA.tx.name = c(circRNA.tx.name, circRNA.as_name)
    }
  }else{ # if not existing the transcript, build circRNA from the exon contig between two ends
    mytx="NONE"
    #get all exons completely inside the boundary
    pick=myexonInfo$EXONSTART >= start & myexonInfo$EXONEND <= end
    myex=myexonInfo[pick,]
    #select the coding area inside the region. 
    #NOTE: this can not always satisfy the canonical splicing conditions (GT-AG), but we do not care about it this momment
    myex=myex[order(myex$EXONSTART),]
    exstart=myex$EXONSTART/1e6
    exend=myex$EXONEND/1e6
    #get clusters of exons
    exCluster=rep(-1,length(exstart))
    for (j in 1:length(exstart)){
      if (exCluster[j]==-1){
        exCluster[j]=j
        pick=(exstart[j]-exstart)*(exstart[j]-exend)<=0 | (exend[j]-exstart)*(exend[j]-exend)<=0
        #assign to the cluster with min index
        x=exCluster[pick]
        x=x[x>0]
        x=min(x)
        exCluster[pick]=x
      }
    }
    #update exCluster
    exCluster_u=exCluster
    repeat{
      exCluster=exCluster[exCluster]
      if (sum(exCluster_u!=exCluster)==0) break()
      exCluster_u=exCluster
    }

    #create exon regions from the exon clusters
    clusterID=sort(unique(exCluster))
    startR=endR=NULL
    for (j in 1:length(clusterID)){
      cID=clusterID[j]
      cEx=myex[exCluster==cID,]
      startR=c(startR,min(cEx$EXONSTART))
      endR=c(endR,max(cEx$EXONEND))
    }

    #sort again startR
    myorder=order(startR)
    startR=startR[myorder]
    endR=endR[myorder]

    # generate pseudo circularRNAs sequences
    #extract sequences and create a new transcript sequence
    chrID=as.character(myex$EXONCHROM)[1]
    x=substring(fasta_genome[chrID],startR,endR)
    txSeq=paste(x,collapse="")
    # new create a pseudo-circularRNA sequence by adding l-1 nucleotides from the end of the last exon to the begin of the first exon
    l_1=readLen -1 -1
    if (l_1 > nchar(txSeq)) l_1=nchar(txSeq) - 1 #if the circRNA length < l_1
    circRNA.seq=paste(substring(txSeq, nchar(txSeq)-l_1),txSeq,sep="")
    circRNA.as_name=paste(mytx,"++",sep="") #do not forget to put the "circRNA", otherwises the format of the output file rawcount.txt of FuSeq is broken  
    #circRNA.fa=circRNA.seq
    #return(c(circRNA.fa,circRNA.as_name))  
    #add to result
    circRNA.tx.fa=c(circRNA.tx.fa,circRNA.seq)
    circRNA.tx.name = c(circRNA.tx.name, circRNA.as_name)
  }
  
  names(circRNA.tx.fa) = circRNA.tx.name
  x = NULL
  circ.name = NULL
  if(length(circRNA.tx.name) > 0){
    x=tapply(circRNA.tx.name, circRNA.tx.fa,c)

    circ.name = sapply(x, function(n){
      z=paste(n, collapse="=")
      paste(circRNA_name, "__", z, " AScircRNA", sep="")
    })
  }

  names(circ.name) = NULL
  return(list("fa" = names(x), "name" = circ.name))
}




#02 Otc 2019
##add get tandem from cDNA for faster getting tandem seq

getTandemSequence_cDNA <-function (myexonInfo, circRNA_name, start, end, num=-1,readLen){
### generate possible tandem sequences for a back-splicing candidate: (start,end) in myexonInfo
  # get sequence
  #output
  tandem.tx.fa = c() #start with a dumpy DNA sequence - serial mode
  tandem.tx.name = c()
  #get all exons completely inside the boundary
  pick1=myexonInfo$EXONSTART == start
  pick2=myexonInfo$EXONEND == end
  tandem_tx=intersect(unique(myexonInfo[pick1,]$TXNAME), unique(myexonInfo[pick2,]$TXNAME))
  #find all transcripts and generate circRNA from those transcripts
  txNum=length(tandem_tx)

 if (txNum > 0){ #Only genertate tandem duplication if existing transcripts containing both EXONSTART and EXONEND    
    #mytx=sample(tandem_tx,1)
    mytx=tandem_tx
    for (atx in mytx){ #do for every tx
      mytxex=myexonInfo[myexonInfo$TXNAME==atx,]
      mytxex=mytxex[order(mytxex$EXONSTART),]
      mytxex.len = mytxex$LENGTH # 
      #covert to transcript location
      mytxex.end.all =cumsum(mytxex$LENGTH)
      mytxex.stat.all = mytxex.end.all - mytxex.len +1
    
      #select the coding area inside the circRNA region
      pick= which(mytxex$EXONSTART >= start & mytxex$EXONEND <= end) #pick by indexes of exons # 
      startR= mytxex.stat.all[pick]
      endR= mytxex.end.all[pick]
          
      if( mytxex$EXONSTRAND[1] == "-"){
      #pick transcript sequence - do converting for - strand
      transcript = convertReverseComplement(as.character(tx.all.fasta[which(tx.all.NAME == atx)]))
      }else{# keep original for + strand
       transcript = as.character(tx.all.fasta[which(tx.all.NAME == atx)])
      }
            
      x=substring(transcript,startR,endR)
      circSeq=paste(x,collapse="") # 
      #pick the left side
      pick.Left=which(mytxex$EXONEND < start) #pick by indexes
      x=""  
      start.Left= mytxex.stat.all[pick.Left]
      end.Left= mytxex.end.all[pick.Left]
      if(length(pick.Left)>0) x = substring(transcript,start.Left,end.Left)
      tandem_left.seq=paste0(x,collapse ="") 
      #pick the right side
      pick.Right=which(mytxex$EXONSTART > end)#pick by indexes
      x=""  
      start.Right= mytxex.stat.all[pick.Right]
      end.Right= mytxex.end.all[pick.Right]
      if(length(pick.Right)>0) x = substring(transcript,start.Right,end.Right)
      tandem_right.seq=paste0(x,collapse ="") # 
      #record tandem break point
      ##paste left side + circSeq and count nchar()
      A.seq = paste0(tandem_left.seq,circSeq)
      #tamdem_break[i] = as.numeric(nchar(A.seq))
      dup_point = as.numeric(nchar(A.seq))+1
      #tandem_normID = my_tandem$normID[1]
      tandem_stat = as.numeric(nchar(tandem_left.seq))+1
      tandem_end = dup_point + as.numeric(nchar(circSeq))
      # finish tandem.seq by pasting again circSeq + tandem_right.seq
      tandem.seq = paste0(A.seq,circSeq,tandem_right.seq) # 
      ##convert back to the correct strand for - strand 
      #if( mytxex$EXONSTRAND[1] == "-") tandem.seq=convertReverseComplement(tandem.seq)
      
      #create tandemRNA.fa and return tadem.tx.fa
      tandem.name=paste(atx,"--",tandem_stat,"--",dup_point,"--",tandem_end, sep="")
      tandem.tx.fa = c(tandem.tx.fa,tandem.seq)
      tandem.tx.name = c(tandem.tx.name, tandem.name)
    }
  }else{ #do nothing now    
    #tandem_normID[i] = "No_transcript"
    #tandem_TXNAME[i] =  "No_transcript"
  }
  names(tandem.tx.fa) = tandem.tx.name
  x = NULL
  tandem.name = NULL
  if(length(tandem.tx.name) > 0){
    x = tapply(tandem.tx.name, tandem.tx.fa, c)
    tandem.name = sapply(x, function(n){
      z=paste(n, collapse="=")
      paste(circRNA_name, "__", z, " tandem", sep="")
    })
  }

  names(tandem.name) = NULL
  return(list("fa" = names(x), "name" = tandem.name))

}



get_circLen <- function(genes.exon.all,circInfo){
  
  genes.exon.all$LENGTH =genes.exon.all$EXONEND - genes.exon.all$EXONSTART +1
  circInfo$normID = paste(circInfo$Chr, circInfo$start_EXONSTART,circInfo$end_EXONEND, sep="__")
  par=foreach(i = 1:nrow(circInfo),.combine=rbind) %dopar% {
    #par=foreach(i = 1:100,.combine=rbind) %dopar% {
    #  if (i %% 500 == 0) cat("\n ",i, " processed.")
      #2) determine the structure of circularRNAs from the list of exons
      mycirR=circInfo[i,]
      mynormID = mycirR$normID
      mygene=as.character(mycirR$GENEID)
      myexonInfo=genes.exon.all[genes.exon.all$GENEID==mygene,] 

      start=mycirR$start_EXONSTART
      end=mycirR$end_EXONEND

      pick1=myexonInfo$EXONSTART == start
      pick2=myexonInfo$EXONEND == end
      mytx=intersect(unique(myexonInfo[pick1,]$TXNAME), unique(myexonInfo[pick2,]$TXNAME))
      #find all transcripts and generate circRNA from those transcripts
      txNum=length(mytx)

      mycircLen = exon_num = start_exon_len = end_exon_len = NULL
      if (txNum > 0){
        for (atx in mytx){ #do for every tx
          mytxex=myexonInfo[myexonInfo$TXNAME==atx,]
          mytxex=mytxex[order(mytxex$EXONSTART),]
          mytxex.len = mytxex$LENGTH
          pick= which(mytxex$EXONSTART >= start & mytxex$EXONEND <= end)
          AS.len = sum(mytxex.len[pick])
          mycircLen = c(mycircLen,AS.len)
          a_exon_num=length(pick)
          start_exon_len=mytxex.len[pick[1]]
          end_exon_len=mytxex.len[pick[a_exon_num]]
          exon_num=c(exon_num,a_exon_num)

        }
      }
      if(txNum == 0){
          #mycircLen=c(0)

          pick=myexonInfo$EXONSTART >= start & myexonInfo$EXONEND <= end
          myex=myexonInfo[pick,]
          #select the coding area inside the region. 
          #NOTE: this can not always satisfy the canonical splicing conditions (GT-AG), but we do not care about it this momment
          myex=myex[order(myex$EXONSTART),]
          exstart=myex$EXONSTART/1e6
          exend=myex$EXONEND/1e6

          if(length(exstart)==0){mycircLen=0}else{
          #get clusters of exons
          exCluster=rep(-1,length(exstart))
          for (j in 1:length(exstart)){
            if (exCluster[j]==-1){
              exCluster[j]=j
              pick=(exstart[j]-exstart)*(exstart[j]-exend)<=0 | (exend[j]-exstart)*(exend[j]-exend)<=0
              #assign to the cluster with min index
              x=exCluster[pick]
              x=x[x>0]
              x=min(x)
              exCluster[pick]=x
            }
          }
          #update exCluster
          exCluster_u=exCluster
          repeat{
            exCluster=exCluster[exCluster]
            if (sum(exCluster_u!=exCluster)==0) break()
            exCluster_u=exCluster
          }     

          #create exon regions from the exon clusters
          clusterID=sort(unique(exCluster))
          startR=endR=NULL
          for (j in 1:length(clusterID)){
            cID=clusterID[j]
            cEx=myex[exCluster==cID,]
            startR=c(startR,min(cEx$EXONSTART))
            endR=c(endR,max(cEx$EXONEND))
          }     

          #sort again startR
          myorder=order(startR)
          startR=startR[myorder]
          endR=endR[myorder]      

          # generate pseudo circularRNAs sequences
          #extract sequences and create a new transcript sequence
          #chrID=as.character(myex$EXONCHROM)[1]
          #x=substring(fasta_genome[chrID],startR,endR)
          #txSeq=paste(x,collapse="")
          exon_len = endR-startR+1
          mycircLen = sum(exon_len)

          exon_num=length(exon_len)
          start_exon_len=exon_len[1]
          end_exon_len=exon_len[exon_num]
          }
        }
    res0 = c(mynormID,min(mycircLen),mean(mycircLen), median(mycircLen),max(mycircLen),txNum,start_exon_len,end_exon_len,max(exon_num)) 
  }

  res = as.data.frame(par,stringsAsFactors = FALSE)
  colnames(res) = c("normID","min.circlen","mean.circlen","median.circlen","max.circlen","txNum","start_exon_len","end_exon_len","max_exon_num")
  res$normID = as.character(res$normID)

  res$min.circlen = as.numeric(res$min.circlen)
  res$mean.circlen = as.numeric(res$mean.circlen)
  res$median.circlen = as.numeric(res$median.circlen)
  res$max.circlen = as.numeric(res$max.circlen)
  res$txNum = as.numeric(res$txNum)
  res$start_exon_len = as.numeric(res$start_exon_len)
  res$end_exon_len = as.numeric(res$end_exon_len)
  res$max_exon_num = as.numeric(res$max_exon_num)


  #out = list(par=par,res=res)
  return(res)
}

process_Sqlite <- function(gtfSqlite){
 anntxdb <- loadDb(gtfSqlite)
 genes.all=genes(anntxdb,single.strand.genes.only=FALSE )
 genes.exon.all= suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID","TXNAME","EXONID","EXONSTART","EXONEND","EXONSTRAND","EXONCHROM"), keytype = "GENEID")))
 #create exon-length # to get seqs from cDNA
 genes.exon.all$LENGTH =genes.exon.all$EXONEND - genes.exon.all$EXONSTART +1
 return(genes.exon.all)
}

process_genome <- function(genomeFastaFile){
  #load genome fasta
  fasta_genome = readDNAStringSet(genomeFastaFile)
  #format the name as 1,2,..X,Y,M..
  chnames=sapply(names(fasta_genome), function(x) unlist(strsplit(x," "))[1])
  #replace MT by M, so it is consistent with refseq annotation
  #chnames[which(chnames=="MT")]="M"
  names(fasta_genome)=chnames
  return(fasta_genome)
}


get_circInfo = function(circRNAs_candidates){
#parse circRNA names to get information
circInfo=strsplit(circRNAs_candidates,"__")
circInfo=do.call(rbind,circInfo)

colnames(circInfo)=c("Chr","start_EXONSTART","end_EXONEND","GENEID","start_EXONID","end_EXONID")
circInfo=as.data.frame(circInfo)

circInfo$start_EXONSTART=as.integer(as.character(circInfo$start_EXONSTART))
circInfo$end_EXONEND=as.integer(as.character(circInfo$end_EXONEND))
return(circInfo)
}

########################################################### for fdr2d

get_fdr2d <- function(in_fn,pre_fn){
  load(in_fn)
  load(pre_fn)

  circRNA_list=PEres$eqc_dat
  nread.obs=circRNA_list$fragment_count
  len.obs=circRNA_list$median.circlen

  set.seed(2019)
  fdr = circfdr(len_obs=len.obs, nread_obs=nread.obs,depdata=Circall_depdata)
  circRNA_list=cbind(circRNA_list,fdr=fdr$ifdr)

  #export to file
  selectColumns=c("normID","fragment_count","fragment_count_norm","median.circlen","fdr")
  circRNA_list=circRNA_list[,selectColumns]
  circInfo=cbind(circInfo,circRNA_list)
  circInfo=circInfo[order(circInfo$fdr,decreasing = FALSE),]
  colnames(circInfo)=c("chr","start","end","geneID","exonID_start","exonID_end","circID","junction_fragment_count", "junction_FPM","median_circlen","fdr")
  circInfo=circInfo[,-c(5,6)]
  return(list(circInfo=circInfo,fdr=fdr))
}



circfdr <- function(len_obs, nread_obs,depdata, readLen=100,
     seed=2019, nrep=100, qlim= 0.999, nr =15, 
     n0=0.1, smooth = 0.9, verb = TRUE, ...){
## len_obs, nread_obs = length and number of supporting fragments of input circRNAs
## depdata = prepared data with column "len" (length), "fragment_count" (number of
##              supporting fragments) and "Type". Column "Type" contain the labels of
##              circRNAs including non_dep and dep
##
## seed = seed for set.seed() 
## nrep = number of random replicates of the null/non-mutated datasets
## qlim = quantile limit of the total reads to avoid outliers
## 
## nr = number of grid points: 2D-area is split to nr x nr
## n0 = proportion of points near the origin where the null is enforced (so fdr=1)
## smooth = smoothing parameter
##
## output: 
## (x,y,fdr.xy) triplet of fdr2d  on the nrxnr 2D-region
## ifdr = individual site-cell fdr estimate: a matrix of mut.sites x ncell
##  

rownames(depdata)=NULL
len_obs=log2(len_obs)
nread_obs=log2(nread_obs)

## null statistics: random sample from presumed not non-depleted circRNAs
#  ## note the seed is used here
if (!is.null(seed)) set.seed(seed) 
n=nrow(depdata)
ran = sample(1:n,n*nrep, replace=TRUE)

permuted_dat=depdata[ran,]
len_all = permuted_dat$len
len_null = permuted_dat$len[permuted_dat$Type != "non_dep"]

nread_all= permuted_dat$fragment_count
nread_null = permuted_dat$fragment_count[permuted_dat$Type != "non_dep"]

x = log2(nread_all)
xnull= log2(nread_null)
y = log2(len_all)
ynull = log2(len_null)

## break points
   xbreaks = MidBreaks(x[x<quantile(x,qlim, na.rm=TRUE) & x>quantile(x,1-qlim, na.rm=TRUE)], nr)
   ybreaks = MidBreaks(y[y<quantile(y,qlim, na.rm=TRUE) & y>quantile(y,1-qlim, na.rm=TRUE)], nr)
   xmid = xbreaks[-1] - diff(xbreaks)/2
   ymid = ybreaks[-1] - diff(ybreaks)/2

## tabulated values
    cut.stat.x = cut(x, xbreaks, include.lowest = TRUE)
    cut.stat.y = cut(y, ybreaks, include.lowest = TRUE)
    cut.null.x = cut(xnull, xbreaks, include.lowest = TRUE)
    cut.null.y = cut(ynull, ybreaks, include.lowest = TRUE)
    tab.stat = table(cut.stat.x, cut.stat.y)
    tab.null = table(cut.null.x, cut.null.y)
 
### need to equalize near the (0,0) point to get fdr=1
   N0 = n0*nr; N0 = ceiling(N0)
   tab.stat[1:N0,1:N0]= tab.null[1:N0,1:N0]= (tab.stat[1:N0,1:N0]+ tab.null[1:N0,1:N0])/2   

## smoothing fdr2d
    df.smooth = (1 - smooth) * (nr - 1)^2
    b = smooth2d.yn(tab.null, tab.stat, df = df.smooth)
    ## note: bfit ~y/n ~ nrep*f0/(nrep*f0 + f)
    ## ==> nrep*f0 = bfit*nrep*f0+ bfit*f)
    ## ==> f0/f = bfit/(1-bfit)/nrep
    f0fz = b$fit/(1 - b$fit)/nrep
    fdr = f0fz
    #fdr = ifelse(f0fz>1,1,f0fz) # bound fdr to 1
    dimnames(fdr) = dimnames(tab.stat)


## indiv fdr values, if needed: cut the fc to reduce computation
    pick = len_obs>= 0
    ifdr = rep(1,length(len_obs))
    ifdr[pick] = approx2d(xmid, ymid, fdr, nread_obs[pick], len_obs[pick])
    ifdr.mat=ifdr

## output
    out = list(x= xmid, y=ymid, fdr.xy=fdr, ifdr=ifdr.mat)
    return(out)

}


## Confidence of a real zero of alt-allele, 
##    accounting for coverage and errors
##
## conf0 = -10*log10(upper-confidence bound of prob)
## based on observing x successes out of n trials
## Example:
#> conf0(0,5); conf0(0,10); conf0(0,30)
#[1] 3.460869 [1] 5.869316 [1] 10.22103
#
conf0 = function(x,n,alpha=0.05){
fn = function(th,x,n, alpha){
    pleft=  pbinom(x,n,th) 
    return(pleft-alpha)
  }
 if (is.na(x)) return(0)
 if (x >= n)  return(0)
 run = uniroot(fn,c(x/n,1),x=x,n=n, alpha=alpha)
 return(-10*log10(run$root))
}



# from OCplus
MidBreaks <- function (x, nr) 
{
    mids = seq(min(x), max(x), length = nr - 1)
    d = mids[2] - mids[1]
    breaks = c(mids - d/2, mids[nr - 1] + d/2)
    breaks
}

## adapted from OCplus simul2d.direct, but with a modification on n==0
## smoothing y/n= proportion of nulls
## NOTE: if n==0 set y/n=1  (all null, so fdr=1)
##
smooth2d.yn<- 
function (y, n, df = 100, err = 0.01, edge.count = 3, mid.start = 6) 
{
    nx = nrow(y)
    ny = ncol(y)
    if (nx < mid.start | edge.count > ny) 
        stop("Grid size too small")
    ## n = ifelse(n <= 1, 1, n)  ## wrong default here!
    edge = 1:edge.count
    mid = mid.start:(nx - mid.start + 1)
    y[mid, edge] = ifelse(n[mid, edge] == 1, 1, y[mid, edge])
    p = y/n
      p = ifelse(n<1, 1, p) ## if n==0 set p=1  (results= null)
    sv2 = df2var(nx, df = df)
    fit = smooth2d.basic(p, sv2 = sv2, err = err)
    list(fit = fit, df = df, sv2 = sv2)
}

## from OCplus
approx2d<- 
function (x, y, z, xout, yout, ...) 
{
    n = length(x)
    m = length(y)
    no = length(xout)
    zout = rep(NA, no)
    for (i in 2:n) {
        for (j in 2:m) {
            ndx = x[i - 1] <= xout & xout <= x[i] & y[j - 1] <= 
                yout & yout <= y[j]
            lower = z[i - 1, j - 1] + (xout[ndx] - x[i - 1]) * 
                (z[i, j - 1] - z[i - 1, j - 1])/(x[i] - x[i - 
                1])
            upper = z[i - 1, j] + (xout[ndx] - x[i - 1]) * (z[i, 
                j] - z[i - 1, j])/(x[i] - x[i - 1])
            zout[ndx] = lower + (yout[ndx] - y[j - 1]) * (upper - 
                lower)/(y[j] - y[j - 1])
        }
    }
    undef = which(is.na(zout))
    if (length(undef > 0)) {
        gg = as.matrix(expand.grid(x, y))
        for (i in undef) {
            d2 = (xout[i] - gg[, 1])^2 + (yout[i] - gg[, 2])
            nn = which.min(d2)
            zout[i] = z[nn]
        }
    }
    zout
}

smooth2d.basic <- function (y, sv2 = 0.01, se2 = 1, err = 0.01) 
{
    W = 1/se2
    nx = nrow(y)
    ny = ncol(y)
    Y = c(y)
    R = rmatrix(nx, ny)
    d = R$d/sv2 + W
    xx = R$x
    ybar = mean(Y)
    B = W * (Y - ybar)
    v = gausei(d, xx, B, sv2, err = err, maxiter = 100)
    est = ybar + v
    matrix(est, ncol = ny)
}

gausei <- function (d, xx, B, sv2, err = 0.01, maxiter = 10) 
{
    ng = length(d)
    xx[xx < 1] = ng + 1
    v = c(B/(d + 1e-06), 0)
    vold = 3 * v
    iter = 1
    while ((sd(vold - v, na.rm = TRUE)/sd(v, na.rm = TRUE) > 
        err) & (iter < maxiter)) {
        vold = v
        vmat = matrix(v[xx], ncol = 4)
        v[1:ng] = (B + c(vmat %*% rep(1, 4))/sv2)/d
        iter = iter + 1
    }
    v[1:ng]
}

rmatrix <- function (nx, ny) 
{
    xymat = matrix(0, nrow = nx, ncol = ny)
    grid = c((row(xymat) - 1) * 1000 + col(xymat) - 1)
    neighbor = c(grid + 1000, grid - 1000, grid + 1, grid - 1)
    loc = match(neighbor, grid, 0)
    loc = matrix(loc, ncol = 4)
    num = c((loc > 0) %*% rep(1, 4))
    list(x = loc, d = num)
}

df2var <- function (nx, df = nx, lower = 1e-04, upper = 1000) 
{
    ff = function(x) dgfree(nx, x) - df
    out = uniroot(ff, lower = lower, upper = upper)
    return(out$root)
}

dgfree <- function (nx, sv2 = 1, ny = nx) 
{
    imp = matrix(0, nrow = nx, ncol = ny)
    imp[nx/2, ny/2] = 1
    imp.resp = smooth2d.basic(imp, sv2 = sv2)
    df = nx * ny * imp.resp[nx/2, ny/2]
    return(df)
}

## OC Curve: true= 0-1 value, pred = probability of true
## extra = extra condition such as VAR_T<0.5
OC = function(true,pred, extra=TRUE, prob=seq(0.01,.99,by=0.01)){
  CUT = SENS= FDR=NULL
  for (pr in prob){
    a = table(true, (pred>pr & extra))
    if (min(dim(a))==1) next
    sens = a[2,2]/(a[2,1]+ a[2,2])
    fdr = a[1,2]/(a[1,2]+a[2,2])
    if ( (fdr-sens)>0 ) next
    SENS = c(SENS, sens)
    FDR = c(FDR,fdr)
    CUT = c(CUT, pr)
  }
  return(list(x=FDR, y=SENS, cut=CUT))
}




