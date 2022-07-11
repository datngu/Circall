
rm(list=ls())
# process arguments and load functions
source("/path/to/R/load_packages.R", print.eval = FALSE)

in_fn=outFn_PE_filtering_Rdata
out_fn=outFn_circRNA_final

#prepare data
load(in_fn)
circRNA_list=PEres$eqc_dat
selectColumns=c("normID","fragment_count","fragment_count_norm", "median.circlen")
circRNA_list=circRNA_list[,selectColumns]
circInfo=cbind(circInfo,circRNA_list)
circInfo=circInfo[order(circInfo$fragment_count,decreasing = TRUE),]
colnames(circInfo)=c("chr","start","end","geneID","exonID_start","exonID_end","circID","junction_fragment_count","junction_FPM","median_circlen")
circInfo=circInfo[,-c(5,6)]
circInfo
#export to file
circRNA_final=circInfo
write.table(circRNA_final,file=out_fn,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
cat("\n Done! Results are save in ",out_fn)

rm(list=ls())
#done! 

PEres$fragInfo_wt$numObservedFragments*10^6