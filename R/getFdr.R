
rm(list=ls())
# process arguments and load functions
source("/path/to/R/load_packages.R", print.eval = FALSE)

in_fn=outFn_PE_filtering_Rdata
out_fn=outFn_circRNA_final
dep_fn=depDataFile

#compute fdr2d
res=get_fdr2d(in_fn,dep_fn)

#export to file
circRNA_final=res$circInfo
write.table(circRNA_final,file=out_fn,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
cat("\n Done! Results are save in ",out_fn)

rm(list=ls())
#done! fdr2d
