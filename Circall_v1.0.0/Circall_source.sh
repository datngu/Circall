#!/bin/bash

#input parameters are:
    #genome -- genome in fasta format
    #gtfSqlite -- genome annotation in Sqlite format
    #txFasta -- transcripts (cDNA) in fasta format
    #txIdx -- quasi-index of txFasta
    #bsjIdx -- quasi-index of BSJ reference fasta file
    #dep -- data contain depleted circRNAs
    #read1 -- input read1
    #read2 -- input read2
    #p -- number of thread
    #tag -- tag name of results
    #td -- generation of tandem sequences
    #c -- clean intermediate data
    #o -- output folder
############################################
syntax="./Circall.sh \n\t-genome [genome in fasta format] \n\t-gtfSqlite [genome annotation in Sqlite format] \n\t-txFasta [transcripts (cDNA) in fasta format] \n\t-txIdx [quasi-index of txFasta] \n\t-bsjIdx [quasi-index of BSJ reference fasta file] \n\t-read1 [read1 fastq.gz file] \n\t-read2 [read2 fastq.gz file] \n\t-thread [number of thread] \n\t-tag [tag name] \n\t-td [TRUE/FALSE value, defaut is FALSE] \n\t-c[TRUE/FALSE value, defaut is TRUE] \n\t-o [output folder]"

#add path
#export LC_ALL=C
#export LD_LIBRARY_PATH=/path/to/linux/lib:$LD_LIBRARY_PATH
#export PATH=/path/to/linux/bin:$PATH

echo -e "
--------------------------------------------------------------------------
                     _                     
            _____   (_)   _____  _____   __        __      __
           / ___/  / /   / ___/ / ___/  /  \      / /     / /
          / /__   / /   / /    / /__   / __ \    / /__   / /__
          \___/  /_/   /_/     \___/  /_/  \_\  /_/__/  /_/__/
                                                                                                                                                                            
---------------------------------------------------------------------------
                     Checking arguments ....
---------------------------------------------------------------------------
"


while [[ $# -gt 1 ]]
do
key="$1"

case $key in

     -genome|--genome) 
     genomeFastaFile=$(readlink -f $2)
     shift
     ;; 

     -gtfSqlite|--annotation) 
     gtfSqlite=$(readlink -f $2)
     shift
     ;; 

     -txFasta|--txFasta) 
     txFastaFile=$(readlink -f $2)
     shift
     ;;

     -txIdx|--txIdx) 
     IndexTranscriptome=$(readlink -f $2)
     shift
     ;;

     -bsjIdx|--bsjIdx) 
     IndexBSJ=$(readlink -f $2)
     shift
     ;;

     -dep|--dep) 
     depDataFile=$(readlink -f $2)
     shift
     ;;


     -read1|--read1) 
     inRead1=$(readlink -f $2)
     shift
     ;;

     -read2|--read2) 
     inRead2=$(readlink -f $2)
     shift
     ;;

     -p|--thread) 
     CPUNUM=$2
     shift
     ;; 

     -tag|--tagname)
     tagStr=$2
     shift
     ;;

     -o|--out)
     outDir=$2
     shift
     ;;

     -c|--clean)
     cleanData=$2
     shift
     ;;

     -td|--tandem) 
     tandem=$2
     shift
     ;; 

     *)

esac
shift
done


if [[ -z "$genomeFastaFile" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no genome sequence file !"
   echo ""
   exit
fi

if [[ -z "$gtfSqlite" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no genome annotation file !"
   echo ""
   exit
fi

if [[ -z "$txFastaFile" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no transcript sequence (cDNA) file !"
   echo ""
   exit
fi

if [[ -z "$IndexTranscriptome" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no index directory of the transcript sequences (cDNA) !"
   echo ""
   exit
fi

if [[ -z "$IndexBSJ" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no index directory of back-splicing-junction (BSJ) reference !"
   echo ""
   exit
fi


if [[ -z "$inRead1" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no input read1 file !"
   echo ""
   exit
fi


if [[ -z "$inRead2" ]]; then
   echo ""
   echo "Usage:"
   echo -e $syntax
   echo ""
   echo "ERROR: no input read2 file !"
   echo ""
   exit
fi



if [[ -z "$depDataFile" ]]; then
   echo "No depeted data file, so will not run fdr2d"
fi


if [[ -z "$CPUNUM" ]]; then
   CPUNUM=4
   echo "No specific setting for thread number, so use the default setting (-p 4)"
fi

if [[ -z "$tagStr" ]]; then
   tagStr=$(echo "Sample")
   echo "No specific setting for tag , so use the default setting (-tag Sample)"
fi

if [[ -z "$tandem" ]]; then
   tandem=$(echo "FASLE")
   echo "No specific setting for --tandem, so use the default setting (-td FASLE)"
fi

if [[ -z "$cleanData" ]]; then
   cleanData=$(echo "TRUE")
   echo "No specific setting for --clean, so use the default setting (-c TRUE)"
fi


if [[ -z "$outDir" ]]; then
   outDir=`pwd`
   echo "No specific setting for output path, so results will be export to current folder"
fi

if [ -d $outDir ]
then
    echo "Checking ... the output folder exists" 
else
    echo "Checking ... the output folder does not exist, so create one"
    mkdir $outDir
fi

in_agr=("genome" "gtfSqlite" "txFasta" "txIdx" "bsjIdx" "dep" "read1" "read2" "thread" "tagname" "tandem" "clean" "output")
used_agrs=("genomeFastaFile" "gtfSqlite" "txFastaFile" "IndexTranscriptome" "IndexBSJ" "depDataFile" "inRead1" "inRead2" "CPUNUM" "tagStr" "tandem" "cleanData" "outDir" "inFileFormat")
inFileFormat=${inRead1##*.}

echo -e "
---------------------------------------------------------------------------
                     Checking arguments successfully!!
---------------------------------------------------------------------------
Your input arguments list:
"

for (( i=0; i< ${#in_agr[@]} ; i++ ));
do
  echo "${in_agr[$i]} was specified as: ${!used_agrs[$i]}"
done

echo ""

# move to the output directory
cd $outDir

root=`pwd`

outDirBS=${root}/${tagStr}_Circall_BS
outDirWT_pairedEnd=${root}/${tagStr}_Circall_wt_UNMAPPED_pairedEnd
outDir_pseudo=${root}/${tagStr}_Circall_pseudo_hitHeader
outDir_pseudo_Idx=${root}/${tagStr}_Circall_pseudo_Idx
Circall_log=${root}/${tagStr}_Circall.log

R_agrs_in="genomeFastaFile=$genomeFastaFile gtfSqlite=$gtfSqlite txFastaFile=$txFastaFile CPUNUM=$CPUNUM"

############################################################################
############################################################################

echo -e "
>>>>>>>>>>>>>>>------------------------------------------------------------
              Step 1: Extracting unmapped reads!!
>>>>>>>>>>>>>>>------------------------------------------------------------
"

# ####### paired end pipeline

Circall_wt -i $IndexTranscriptome -1 <(gunzip -c $inRead1) -2 <(gunzip -c $inRead2) -o $outDirWT_pairedEnd -p $CPUNUM >>$Circall_log 2>&1

unmappedFilesMerged1="$outDirWT_pairedEnd/Unmapped_readM_1.fa"
unmappedFilesMerged2="$outDirWT_pairedEnd/Unmapped_readM_2.fa"

# Merge the outputs from multiple cores
cat $outDirWT_pairedEnd/UN_read1_* >  $unmappedFilesMerged1
rm $outDirWT_pairedEnd/UN_read1_*

cat $outDirWT_pairedEnd/UN_read2_* >  $unmappedFilesMerged2
rm $outDirWT_pairedEnd/UN_read2_*


echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>---------------------------------------------
               Step 2: Mapping to BSJ reference !!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>---------------------------------------------
"
Circall_bsj -i $IndexBSJ -1 $unmappedFilesMerged1 -2 $unmappedFilesMerged2 -o $outDirBS -p $CPUNUM >>$Circall_log 2>&1

#remove unmapped reads to release storage
rm $unmappedFilesMerged1 $unmappedFilesMerged2

#process BSJ reads
BSJ_read_merged_1="$outDirBS/BSJ_read_merged_1.fa"
BSJ_read_merged_2="$outDirBS/BSJ_read_merged_2.fa"

# Merge the outputs from multiple cores
cat $outDirBS/BSJ_read1_* > $BSJ_read_merged_1
rm $outDirBS/BSJ_read1_*

cat $outDirBS/BSJ_read2_* > $BSJ_read_merged_2
rm $outDirBS/BSJ_read2_*

echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
              Step 3: Filtering for single-end reads!!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
" 

outFn_SE_filtering_Rdata=${root}/${tagStr}_Circall_SE_filter_output.RData
# do the SE filter
Rscript /path/to/R/doSingleEndFiltering.R $R_agrs_in outDirBS=$outDirBS outFn_SE_filtering_Rdata=$outFn_SE_filtering_Rdata >>$Circall_log 2>&1

echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>---------------
              Step 4: Generating pseudo transcripts !!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>---------------
" 

#extract putative circRNA sequences base on SE filtering results
outFn_getPseudoSeq_Rdata=${root}/${tagStr}_Circall_circRNAinfo_pseudoSeq.RData
outFn_getPseudoSeq_fasta=${root}/${tagStr}_Circall_pseudoSeq.fa


Rscript /path/to/R/getPseudoCircRNAsequence.R $R_agrs_in outDirBS=$outDirBS outFn_SE_filtering_Rdata=$outFn_SE_filtering_Rdata outFn_getPseudoSeq_Rdata=$outFn_getPseudoSeq_Rdata outFn_getPseudoSeq_fasta=$outFn_getPseudoSeq_fasta tandem=$tandem >>$Circall_log 2>&1


echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              Step 5: Mapping paired-end reads!!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
" 

TxIndexer -t $outFn_getPseudoSeq_fasta -p $CPUNUM -o $outDir_pseudo_Idx -f >>$Circall_log 2>&1

##map paired-end and extract the header of the hits
Circall_pseudo -i $outDir_pseudo_Idx -1 $BSJ_read_merged_1 -2 $BSJ_read_merged_2 -o $outDir_pseudo -p $CPUNUM >>$Circall_log 2>&1
##merge
cat $outDir_pseudo/*.rh > $outDir_pseudo/merged_hitheader.txt
rm $outDir_pseudo/*.rh
cat $outDir_pseudo/mapInfo_*.txt > $outDir_pseudo/merged_mapInfo.txt
rm $outDir_pseudo/mapInfo_*.txt


###PE filter and do filter for circRNA.AS fasta file#
outFn_PE_filtering_Rdata=${root}/${tagStr}_Circall_PE_filter_output.RData

Rscript /path/to/R/doPairEndFiltering.R $R_agrs_in outFn_SE_filtering_Rdata=$outFn_SE_filtering_Rdata outDirWT=$outDirWT_pairedEnd outDir_pseudo=$outDir_pseudo outFn_PE_filtering_Rdata=$outFn_PE_filtering_Rdata >>$Circall_log 2>&1


outFn_circRNA_final=${root}/${tagStr}_Circall_final.txt
if [[ -z "$depDataFile" ]]; then
echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                Do not compute fdr2d, just export to file!!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
"
Rscript /path/to/R/export.R outFn_PE_filtering_Rdata=$outFn_PE_filtering_Rdata outFn_circRNA_final=$outFn_circRNA_final >>$Circall_log 2>&1

else

echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                Compute fdr2d and export to file!!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
"
Rscript /path/to/R/getFdr.R outFn_PE_filtering_Rdata=$outFn_PE_filtering_Rdata depDataFile=$depDataFile outFn_circRNA_final=$outFn_circRNA_final >>$Circall_log 2>&1
fi


if [[ "$cleanData" == "TRUE" ]]; then
echo -e "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                Clean all intermediate data !!
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
" 
rm -r $outDirWT_pairedEnd
rm -r $outDirBS
rm -r $outFn_getPseudoSeq_fasta
rm -r $outDir_pseudo_Idx
rm -r $outDir_pseudo
fi

echo -e "

--------------------------------------------------------------------------
                     _                     
            _____   (_)   _____  _____   __        __      __
           / ___/  / /   / ___/ / ___/  /  \      / /     / /
          / /__   / /   / /    / /__   / __ \    / /__   / /__
          \___/  /_/   /_/     \___/  /_/  \_\  /_/__/  /_/__/

-------------------------------------------------------------------------
                            DONE!
---------------------------------------------------------------------------
" 
