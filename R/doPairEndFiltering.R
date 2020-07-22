# This module used to get pseudo circRNA sequences from results of SE filtering

# Input: outFn_SE_filtering_Rdata, path/to/outDir_pseudo, path/to/outDirWT
# Output: outFn_PE_filtering_Rdata that contains CircRNA list and needed information for fdr2d

# Syntax doPairEndFiltering.R outDirWT=path/to/outDirWT outDir_pseudo=path/to/outDir_pseudo outFn_SE_filtering=path/to/SE_filtering.Rdata outFn_PE_filtering=path/to/PE_filtering.Rdata CPUNUM=cpu_number



rm(list = ls())

source("/path/to/R/load_packages.R", print.eval = FALSE)

# load results from outDir_pseudo
path_to_header =  paste0(outDir_pseudo, "/merged_hitheader.txt" )
path_to_mapInfo =  paste0(outDir_pseudo, "/merged_mapInfo.txt" )
PEread = load_read2_2_file(path_to_header,path_to_mapInfo)

# take fragmentInfo from outDirWT, moved from SE_filtering script
fragInfo_wt=read.csv(paste0(outDirWT,"/fragmentInfo.txt"),sep="\t")

# load results from SE_filtering
SE_filtering_res=load_SE_filtering(outFn_SE_filtering_Rdata)

# do Pair-end filtering
PEres = PE_filter(PEread = PEread, resPhase1 = SE_filtering_res, fragInfo_wt = fragInfo_wt, read_mapProp = 0, read_anchorLen = 10, eqc_txMinProp = 1, min_fragment_count = 1)

circInfo = get_circInfo(as.character(PEres$circRNAs_candidates))

## This part for Xtest - Skip it now!
# path_to_Info =paste(tagStr,"_Circall_circRNAinfo_pseudoSeq.RData",sep="")
# load(path_to_Info) #load info of tx.fa created by previous step

# #filter circRNA.AS
# tx.name = sapply(as.character(names(pseudoSeq.fa)),function(x){
#   y=gregexpr("__",x)
#   #s=y[[1]][3]+2
#   e=y[[1]][6]-1
#   return(substring(x,1,e))
# })
# pseudo_circRNA.AS.tx = pseudoSeq.fa[tx.name %in% PEres$circRNAs_candidates,]

# save(pseudo_circRNA.AS.tx, circInfo, PEres, file = outFn_PE_filtering_Rdata)


save(circInfo, PEres, file = outFn_PE_filtering_Rdata)

##### Done - PE filtering

rm(list=ls())