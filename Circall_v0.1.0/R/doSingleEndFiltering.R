# This module is used to filter based on single end mapping information
# Filtering criteria:
# 1. read_anchorLen=10, minumum read anchor length.
# 2. eqc_txMinProp= 0.99, minumum propotation of read length is mapped to BSJ reference.
# 3. eqc_txMinProp = 1, each read is mapped unique to only 1 equevalent class.
# 4. eqc_minCount=1, minumun number of reads mapped to a equevalent class.
# 5. hard_hitlen=50, for short circRNA, the mimimum mapped length is set with 50.

# Input: a dir path results of Circall_wt, gtfSqlite annotation
# Output: an RData file that contains raw list of circRNA and downtream needed information.
# Syntax doSingleEndFiltering.R outDirBS=path/to/outDirBS outFn_SE_filtering=path/to/SE_filtering.Rdata gtfSqlite=path/to/gtf.sqlite CPUNUM=cpu_number


###################################
rm(list=ls())
# process arguments and load functions
source("/path/to/R/load_packages.R", print.eval = FALSE)
#do filtering
SE_filtering=SE_filtering(inPath = outDirBS, gtfSqlite = gtfSqlite, outFn = outFn_SE_filtering_Rdata)
rm(list=ls())
#done! SE_filtering
