/*
Date: 06 Mar 2020
This tool is used to extract BSJ candidates from RNA-seq data.
Note:This implementation is adapted from SailfishQuantify.cpp of Sailfish for our purpose. 
*/

#include <algorithm>
#include <cstdio>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <thread>
#include <string>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

// TBB include
#include "tbb/atomic.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

//#include "BiasIndex.hpp"
#include "VersionChecker.hpp"
#include "SailfishConfig.hpp"
#include "SailfishUtils.hpp"
#include "SailfishIndex.hpp"
#include "TranscriptGeneMap.hpp"
//#include "CollapsedEMOptimizer.hpp"
#include "ExportFeq.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "ReadLibrary.hpp"
#include "RapMapUtils.hpp"
#include "HitManager.hpp"
#include "SASearcher.hpp"
//#include "SACollector.hpp"
#include "SACollectorCircall.hpp"
#include "EmpiricalDistribution.hpp"
#include "TextBootstrapWriter.hpp"
#include "GZipWriter.hpp"
//#include "HDF5Writer.hpp"

#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"




#include "spdlog/spdlog.h"

//S_AYUSH_CODE
#include "ReadKmerDist.hpp"
//T_AYUSH_CODE

/****** QUASI MAPPING DECLARATIONS *********/
using MateStatus = rapmap::utils::MateStatus;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
/****** QUASI MAPPING DECLARATIONS  *******/

/****** Parser aliases ***/
//using paired_parser = pair_sequence_parser<std::vector<std::ifstream*>::iterator>;
using paired_parser = pair_sequence_parser<char**>;//std::vector<std::ifstream*>::iterator>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
/****** Parser aliases ***/


// using FragLengthCountMap = std::unordered_map<uint32_t, uint64_t>;
using FragLengthCountMap = std::vector<tbb::atomic<uint32_t>>;

using std::string;

constexpr uint32_t readGroupSize{1000};

/**
 * Compute and return the mean fragment length ---
 * rounded down to the nearest integer --- of the fragment
 * length distribution.
 */
int32_t getMeanFragLen(const FragLengthCountMap& flMap) {
    double totalCount{0.0};
    double totalLength{0.0};
    for (size_t i = 0; i < flMap.size(); ++i) {
        auto c = flMap[i];
        totalLength += i * c;
        totalCount += c;
    }
    double ret{200.0};
    if (totalCount <= 0.0) {
        std::cerr << "Saw no fragments; can't compute mean fragment length.\n";
        std::cerr << "This appears to be a bug. Please report it on GitHub.\n";
        return ret;
    }
    if (totalLength > totalCount) {
        ret = (totalLength / totalCount);
    }
    return static_cast<uint32_t>(ret);
}




/**
 * For paired-end reads:
 * Do the main work of mapping the reads and building
 * the equivalence classes.
 */
template <typename IndexT>
void processReadsQuasi(paired_parser* parser,
               IndexT* sidx,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               SailfishOpts& sfOpts,
               FragLengthCountMap& flMap,
               std::atomic<int32_t>& remainingFLOps,
             std::mutex& iomutex) {

  uint32_t maxFragLen = sfOpts.maxFragLen;
  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};
  int32_t meanFragLen{-1};

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};

  bool tooManyHits{false};
  size_t maxNumHits{sfOpts.maxReadOccs};
  size_t readLen{0};

  auto& jointLog = sfOpts.jointLog;
  auto& numObservedFragments = readExp.numObservedFragmentsAtomic();
  auto& validHits = readExp.numMappedFragmentsAtomic();
  auto& totalHits = readExp.numFragHitsAtomic();
  auto& upperBoundHits = readExp.upperBoundHitsAtomic();
  auto& eqBuilder = readExp.equivalenceClassBuilder();

  auto& transcripts = readExp.transcripts();

  auto& readBias = readExp.readBias();
  auto& observedGC = readExp.observedGC();
  bool estimateGCBias = sfOpts.gcBiasCorrect;
  bool strictIntersect = sfOpts.strictIntersect;
  bool discardOrphans = !sfOpts.allowOrphans;

//  SACollector<IndexT> hitCollector(sidx);
  SACollectorCircall<IndexT> hitCollector(sidx);
  
  SASearcher<IndexT> saSearcher(sidx);
  rapmap::utils::HitCounters hctr;

  std::vector<QuasiAlignment> leftHits;
  std::vector<QuasiAlignment> rightHits;
  std::vector<QuasiAlignment> jointHits;


  //Circall
  //template <typename RapMapIndexT>
  using OffsetT = typename IndexT::IndexType;
  using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;
  std::vector<SAIntervalHit> leftFwdSAInts;
  std::vector<SAIntervalHit> leftRcSAInts;
  std::vector<SAIntervalHit> rightFwdSAInts;
  std::vector<SAIntervalHit> rightRcSAInts;

  auto& SA = sidx->SA;
  auto& txpStarts = sidx->txpOffsets;  

  std::vector<uint32_t> txpIDsAll;
  std::vector<double> auxProbsAll;

  std::vector<uint32_t> txpIDsCompat;
  std::vector<double> auxProbsCompat;

  // *Completely* ignore strandedness information
  bool ignoreCompat = sfOpts.ignoreLibCompat;
  // Don't *strictly* enforce compatibility --- if
  // the only hits are incompatible with the library
  // type then allow them.
  bool enforceCompat = sfOpts.enforceLibCompat;
  // True when we have compatible hits, false otherwise
  bool haveCompat{false};
  auto expectedLibType = rl.format();

  bool canDovetail = sfOpts.allowDovetail;

  bool mappedFrag{false};
  std::unique_ptr<EmpiricalDistribution> empDist{nullptr};

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, sfOpts.maxReadOccs);
  //Circall
  auto kmerLen = rapmap::utils::my_mer::k();
  //generate a unique ID for as single thread
  string ftag=std::to_string(std::rand()); 
  string UN_fafn1 = "UN_read1_" + ftag +".fa";
  boost::filesystem::path UN_readfafn1 = sfOpts.outputDirectory  / UN_fafn1;
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> UN_faOutput1(std::fopen(UN_readfafn1.c_str(), "a"), std::fclose);
//  string UN_fafn2 = "UN_read2_" + ftag +".fa";
//  boost::filesystem::path UN_readfafn2 = sfOpts.outputDirectory  / UN_fafn2;
//  std::unique_ptr<std::FILE, int (*)(std::FILE *)> UN_faOutput2(std::fopen(UN_readfafn2.c_str(), "a"), std::fclose);

  string BSJ_fafn1 = "BSJ_read1_" + ftag +".fa";
  boost::filesystem::path BSJ_readfafn1 = sfOpts.outputDirectory  / BSJ_fafn1;
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> BSJ_faOutput1(std::fopen(BSJ_readfafn1.c_str(), "a"), std::fclose);
  string BSJ_fafn2 = "BSJ_read2_" + ftag +".fa";
  boost::filesystem::path BSJ_readfafn2 = sfOpts.outputDirectory  / BSJ_fafn2;
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> BSJ_faOutput2(std::fopen(BSJ_readfafn2.c_str(), "a"), std::fclose);


  while(true) {
    typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
        readLen = j->data[i].first.seq.length();
        tooManyHits = false;
        jointHits.clear();
        leftHits.clear();
        rightHits.clear();
        txpIDsAll.clear();
        auxProbsAll.clear();
        txpIDsCompat.clear();
        auxProbsCompat.clear();
        haveCompat = false;
        mappedFrag = false;

        //Circall
        leftRcSAInts.clear();  
        leftFwdSAInts.clear();
//        rightFwdSAInts.clear();
//        rightRcSAInts.clear();
        
/*
        bool lh = hitCollector(j->data[i].first.seq,
                               leftHits, saSearcher,
                               MateStatus::PAIRED_END_LEFT,
                 true // strict check
                 );

        bool rh = hitCollector(j->data[i].second.seq,
                               rightHits, saSearcher,
                               MateStatus::PAIRED_END_RIGHT,
                 true // strict check
                 );
*/
/*
                //export information  of a read pair
        bool lh = hitCollector(j->data[i].first.seq,
                               leftHits, 
                               leftFwdSAInts, leftRcSAInts, 
                               saSearcher,
                               MateStatus::PAIRED_END_LEFT,
                 true // strict check
                 );

        bool rh = hitCollector(j->data[i].second.seq,
                               rightHits, 
                               rightFwdSAInts, rightRcSAInts, 
                               saSearcher,
                               MateStatus::PAIRED_END_RIGHT,
                 true // strict check
                 );

        if (strictIntersect) {
          rapmap::utils::mergeLeftRightHits(
              leftHits, rightHits, jointHits,
              readLen, maxNumHits, tooManyHits, hctr);
        } else {
          rapmap::utils::mergeLeftRightHitsFuzzy(
              lh, rh,
              leftHits, rightHits, jointHits,
              readLen, maxNumHits, tooManyHits, hctr);
        }
        */
          leftRcSAInts.clear();  
          leftFwdSAInts.clear();

/*
          bool lh = hitCollector(j->data[i].first.seq,
                  jointHits, 
                  leftRcSAInts, leftFwdSAInts,
                  saSearcher,
                  MateStatus::SINGLE_END);
*/
          bool lh = hitCollector(j->data[i].first.seq,
                       jointHits, 
                       leftFwdSAInts, leftRcSAInts, 
                       saSearcher,
//                       MateStatus::PAIRED_END_LEFT,
                       MateStatus::SINGLE_END,
         true // strict check
         );

        //Circall
        if (readExp.getReadLength() == 0) readExp.setReadLength(readLen);
        std::string header1=">"+j->data[i].first.header;
        std::string header2=">"+j->data[i].second.header;

//////////////////////////////

       /*----------------------------------------------------------------*/
        // There are 02 kmer lists. Need to check whether their k-mers contain back-splicing junction or not 
        bool isSplit{false};
        if (jointHits.size()>0){ 
          // Case 1: leftFwdValid that exists hits from forward direction: Detect all cases that hitPos < 100 and hitPos+len>100
          //Check SA intervals
//            std::cout << " leftFwdSAInts, kmer number " << leftFwdSAInts.size()<<std::endl; 
              for (auto& saIntervalHit : leftFwdSAInts) {
//                std::cout << " begin " << saIntervalHit.begin<< " end " << saIntervalHit.end << " len "<< saIntervalHit.len<< " queryPos "<<saIntervalHit.queryPos <<std::endl; 
              //auto& saIntervalHit = leftFwdSAInts.front();
              for (OffsetT sa_i = saIntervalHit.begin; sa_i != saIntervalHit.end; ++sa_i) {
                  auto globalPos = SA[sa_i];
                  auto txpID = sidx->transcriptAtPosition(globalPos);
                  // the offset into this transcript
                  auto pos = globalPos - txpStarts[txpID];
                  int32_t hitPos = pos;
                  //int32_t junctPos = 99;//transcripts[txpID].RefLength/2 = 198/2
                  int32_t junctPos = floor(transcripts[txpID].RefLength/2);
//                  std::cout <<std::endl << " txpID " <<txpID << " txName " <<transcripts[txpID].RefName << " txLen " << transcripts[txpID].RefLength << " hitPos "<< std::to_string(hitPos) ; 
                  //If change something here, do not forget change in case 2                  
                  if ((hitPos < (junctPos-5)) and ((hitPos+saIntervalHit.len) > (junctPos+5))){
//                  if ((hitPos < 89) and ((hitPos + std::stoi(std::to_string(saIntervalHit.len))) > 109)){
                    //FW(0)/RC(1)
                    std::string myJunction=header1+"\t" + "0" +"\t"+std::to_string(saIntervalHit.queryPos)+"\t"+std::to_string(saIntervalHit.len) +"\t"+std::to_string(hitPos)+ "\t"+transcripts[txpID].RefName + "\t"+std::to_string(transcripts[txpID].RefLength);
                    fmt::print(UN_faOutput1.get(), myJunction);
                    fmt::print(UN_faOutput1.get(), "\n");
                    isSplit=true;
                    //std::cout <<std::endl  << myJunction  << std::endl;

                  }
              }
            }

          // Case 2: leftRcSAInts that exists hits from forward direction: Detect all cases that hitPos < 100 and hitPos+len>100
          //Check SA intervals
//            std::cout << " leftRcSAInts, kmer number " << leftRcSAInts.size()<<std::endl; 
              for (auto& saIntervalHit : leftRcSAInts) {
//                std::cout << " begin " << saIntervalHit.begin<< " end " << saIntervalHit.end << " len "<< saIntervalHit.len<< " queryPos "<<saIntervalHit.queryPos <<std::endl; 
              //auto& saIntervalHit = leftRcSAInts.front();
              for (OffsetT sa_i = saIntervalHit.begin; sa_i != saIntervalHit.end; ++sa_i) {
                  auto globalPos = SA[sa_i];
                  auto txpID = sidx->transcriptAtPosition(globalPos);
                  // the offset into this transcript
                  auto pos = globalPos - txpStarts[txpID];
                  int32_t hitPos = pos;                  
                  //int32_t junctPos = 99;//transcripts[txpID].RefLength/2 = 198/2
                  int32_t junctPos = floor(transcripts[txpID].RefLength/2);
//                  std::cout <<std::endl << " txpID " <<txpID << " txName " <<transcripts[txpID].RefName << " txLen " << transcripts[txpID].RefLength << " hitPos "<< std::to_string(hitPos) ; 
                  if ((hitPos < (junctPos-5)) and ((hitPos+saIntervalHit.len) > (junctPos+5))){
                    //FW(0)/RC(1)
                    std::string myJunction=header1+"\t"+ "1" +"\t" + std::to_string(saIntervalHit.queryPos)+"\t"+std::to_string(saIntervalHit.len) +"\t"+std::to_string(hitPos) +"\t"+transcripts[txpID].RefName + "\t"+std::to_string(transcripts[txpID].RefLength);
                    fmt::print(UN_faOutput1.get(), myJunction);
                    fmt::print(UN_faOutput1.get(), "\n");
                    isSplit=true;
                    //std::cout <<std::endl  << myJunction  << std::endl;

                  }
              }
            }


          }
/////////////////////////////




        upperBoundHits += (jointHits.size() > 0);

        if (jointHits.size() > sfOpts.maxReadOccs ) { jointHits.clear(); }


            //Circall: just update equivalence classes if there are split reads
            //if (jointHits.size() > 0) {
            if (jointHits.size() > 0 and isSplit) {
            // export if there exists split reads
            fmt::print(BSJ_faOutput1.get(), header1);
            fmt::print(BSJ_faOutput1.get(), "\n");
            fmt::print(BSJ_faOutput1.get(), j->data[i].first.seq);
            fmt::print(BSJ_faOutput1.get(), "\n");
            fmt::print(BSJ_faOutput2.get(), header2);
            fmt::print(BSJ_faOutput2.get(), "\n");
            fmt::print(BSJ_faOutput2.get(), j->data[i].second.seq);
            fmt::print(BSJ_faOutput2.get(), "\n");

       
                int32_t fwAll = 0;
                int32_t fwCompat = 0;
                int32_t rcAll = 0;
                int32_t rcCompat = 0;

                double auxSumAll = 0.0;
                double auxSumCompat = 0.0;

                bool needBiasSample = sfOpts.biasCorrect;

                for (auto& h : jointHits) {
                    auto transcriptID = h.transcriptID();
                    auto& txp = transcripts[transcriptID];

                    int32_t pos = static_cast<int32_t>(h.pos);
                    auto dir = sailfish::utils::boolToDirection(h.fwd);

                    // Note: sidx is a pointer to type IndexT, not RapMapSAIndex!

                    // If bias correction is turned on, and we haven't sampled a mapping
                    // for this read yet, and we haven't collected the required number of
                    // samples overall.
                    if(needBiasSample and sfOpts.numBiasSamples > 0){
                        // the "start" position is the leftmost position if
                        // we hit the forward strand, and the leftmost
                        // position + the read length if we hit the reverse complement
                        int32_t startPos = h.fwd ? pos : pos + h.readLen;

                        if (startPos > 0 and startPos < txp.RefLength) {
   
                            const char* txpStart = txp.Sequence();
                            const char* readStart = txpStart + startPos; // is this correct?
                            const char* txpEnd = txpStart+ txp.RefLength;
                            bool success = readBias.update(txpStart, readStart, txpEnd, dir);
                            if (success) {
                                sfOpts.numBiasSamples -= 1;
                                needBiasSample = false;
                            }
                        }
                    }

                    // True if the read is compatible with the
                    // expected library type; false otherwise.
                    bool compat = ignoreCompat;
                    if (!compat) {
                        compat = sailfish::utils::compatibleHit(
                                expectedLibType, pos,
                                h.fwd, h.mateStatus);
                    }

                    if (compat) {
                        haveCompat = true;
                        txpIDsCompat.push_back(transcriptID);
                        auxProbsCompat.push_back(1.0);
                        auxSumCompat += 1.0;
                        if (h.fwd) { fwCompat++; } else { rcCompat++; }
                    }
                    if (!haveCompat and !enforceCompat) {
                        txpIDsAll.push_back(transcriptID);
                        auxProbsAll.push_back(1.0);
                        auxSumAll += 1.0;
                        if (h.fwd) { fwAll++; } else { rcAll++; }
                    }
        }

                // If we have compatible hits, only use those
                if (haveCompat) {
                    if (txpIDsCompat.size() > 0) {
                        mappedFrag = true;
                        TranscriptGroup tg(txpIDsCompat);
                        eqBuilder.addGroup(std::move(tg), auxProbsCompat);
                        readExp.addNumFwd(fwCompat);
                        readExp.addNumRC(rcCompat);
                    }
                } else {
                    if (txpIDsAll.size() > 0) {
                        // Otherwise, consider all hits.
                        mappedFrag = true;
                        TranscriptGroup tg(txpIDsAll);
                        eqBuilder.addGroup(std::move(tg), auxProbsAll);
                        readExp.addNumFwd(fwAll);
                        readExp.addNumRC(rcAll);
                    }
                }
            }

            validHits += (mappedFrag) ? 1 : 0;
            totalHits += jointHits.size();
            locRead++;
            ++numObservedFragments;
            if (numObservedFragments % 500000 == 0) {
                iomutex.lock();
                fmt::print(stderr, "\033[A\r\rprocessed {} fragments\n", numObservedFragments);
                fmt::print(stderr, "hits: {}, hits per frag (may not be concordant):  {}",
                        totalHits,
                        totalHits / static_cast<float>(prevObservedFrags));

                iomutex.unlock();
            }

    } // end for i < j->nb_filled
    prevObservedFrags = numObservedFragments;
  }
}


/**
 * For single-end reads:
 * Map the reads and accumulate equivalence class counts.
 **/
template <typename IndexT>
void processReadsQuasi(single_parser* parser,
        IndexT* sidx,
        ReadExperiment& readExp,
        ReadLibrary& rl,
        SailfishOpts& sfOpts,
        std::mutex& iomutex) {

    uint64_t prevObservedFrags{1};

    size_t locRead{0};
    uint64_t localUpperBoundHits{0};
    //S_AYUSH_CODE
    auto& readBias = readExp.readBias();
    const char* txomeStr = sidx->seq.c_str();
    //T_AYUSH_CODE

    bool tooManyHits{false};
    size_t readLen{0};
    size_t maxNumHits{sfOpts.maxReadOccs};

    auto& numObservedFragments = readExp.numObservedFragmentsAtomic();
    auto& validHits = readExp.numMappedFragmentsAtomic();
    auto& totalHits = readExp.numFragHitsAtomic();
    auto& upperBoundHits = readExp.upperBoundHitsAtomic();
    auto& eqBuilder = readExp.equivalenceClassBuilder();
    auto& transcripts = readExp.transcripts();

    //auto sidx = readExp.getIndex();
    SACollectorCircall<IndexT> hitCollector(sidx);
    SASearcher<IndexT> saSearcher(sidx);
    rapmap::utils::HitCounters hctr;
    std::vector<QuasiAlignment> jointHits;

    // *Completely* ignore strandedness information
    bool ignoreCompat = sfOpts.ignoreLibCompat;
    // Don't *strictly* enforce compatibility --- if
    // the only hits are incompatible with the library
    // type then allow them.
    bool enforceCompat = sfOpts.enforceLibCompat;
    // True when we have compatible hits, false otherwise
    bool haveCompat{false};
    auto expectedLibType = rl.format();

    bool mappedFrag{false};

    std::vector<uint32_t> txpIDsAll;
    std::vector<double> auxProbsAll;

    std::vector<uint32_t> txpIDsCompat;
    std::vector<double> auxProbsCompat;

    //Circall
    using OffsetT = typename IndexT::IndexType;
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;
    std::vector<SAIntervalHit> leftFwdSAInts;
    std::vector<SAIntervalHit> leftRcSAInts;
  //  std::vector<SAIntervalHit> rightFwdSAInts;
  //  std::vector<SAIntervalHit> rightRcSAInts;
    auto& SA = sidx->SA;
    auto& txpStarts = sidx->txpOffsets;  
  //generate a unique ID for as single thread
  string ftag=std::to_string(std::rand());
  string UN_fafn1 = "UN_read_" + ftag +".fa";
  boost::filesystem::path UN_readfafn1 = sfOpts.outputDirectory  / UN_fafn1;
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> UN_faOutput1(std::fopen(UN_readfafn1.c_str(), "a"), std::fclose);

    while(true) {
        typename single_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;        // If got nothing, quit

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
            readLen = j->data[i].seq.length();
            tooManyHits = false;
            localUpperBoundHits = 0;
            jointHits.clear();
            txpIDsAll.clear();
            auxProbsAll.clear();
            txpIDsCompat.clear();
            auxProbsCompat.clear();
            haveCompat = false;
            mappedFrag = false;

            //Circall
            leftRcSAInts.clear();  
            leftFwdSAInts.clear();


            bool lh = hitCollector(j->data[i].seq,
                    jointHits, 
                    //leftRcSAInts, leftFwdSAInts,
                    leftFwdSAInts,leftRcSAInts,
                    saSearcher,
                    MateStatus::SINGLE_END);
            //Circall
            if (readExp.getReadLength() == 0) readExp.setReadLength(readLen);
            std::string header1=">"+j->data[i].header;




       /*----------------------------------------------------------------*/
        // There are 02 kmer lists. Need to check whether their k-mers contain back-splicing junction or not 
        bool isSplit{false};
        if (jointHits.size()>0){ 
          // Case 1: leftFwdValid that exists hits from forward direction: Detect all cases that hitPos < 100 and hitPos+len>100
          //Check SA intervals
//            std::cout << " leftFwdSAInts, kmer number " << leftFwdSAInts.size()<<std::endl; 
              for (auto& saIntervalHit : leftFwdSAInts) {
//                std::cout << " begin " << saIntervalHit.begin<< " end " << saIntervalHit.end << " len "<< saIntervalHit.len<< " queryPos "<<saIntervalHit.queryPos <<std::endl; 
              //auto& saIntervalHit = leftFwdSAInts.front();
              for (OffsetT sa_i = saIntervalHit.begin; sa_i != saIntervalHit.end; ++sa_i) {
                  auto globalPos = SA[sa_i];
                  auto txpID = sidx->transcriptAtPosition(globalPos);
                  // the offset into this transcript
                  auto pos = globalPos - txpStarts[txpID];
                  int32_t hitPos = pos;
                  //int32_t junctPos = 99;//transcripts[txpID].RefLength/2 = 198/2
                  int32_t junctPos = floor(transcripts[txpID].RefLength/2);
//                  std::cout <<std::endl << " txpID " <<txpID << " txName " <<transcripts[txpID].RefName << " txLen " << transcripts[txpID].RefLength << " hitPos "<< std::to_string(hitPos) ; 
                  //If change something here, do not forget change in case 2                  
                  if ((hitPos < (junctPos-5)) and ((hitPos+saIntervalHit.len) > (junctPos+5))){
//                  if ((hitPos < 89) and ((hitPos + std::stoi(std::to_string(saIntervalHit.len))) > 109)){
                    //FW(0)/RC(1)
                    std::string myJunction=header1+"\t" + "0" +"\t"+std::to_string(saIntervalHit.queryPos)+"\t"+std::to_string(saIntervalHit.len) +"\t"+std::to_string(hitPos)+ "\t"+transcripts[txpID].RefName + "\t"+std::to_string(transcripts[txpID].RefLength);
                    fmt::print(UN_faOutput1.get(), myJunction);
                    fmt::print(UN_faOutput1.get(), "\n");
                    isSplit=true;
                    //std::cout <<std::endl  << myJunction  << std::endl;

                  }
              }
            }

          // Case 2: leftRcSAInts that exists hits from forward direction: Detect all cases that hitPos < 100 and hitPos+len>100
          //Check SA intervals
//            std::cout << " leftRcSAInts, kmer number " << leftRcSAInts.size()<<std::endl; 
              for (auto& saIntervalHit : leftRcSAInts) {
//                std::cout << " begin " << saIntervalHit.begin<< " end " << saIntervalHit.end << " len "<< saIntervalHit.len<< " queryPos "<<saIntervalHit.queryPos <<std::endl; 
              //auto& saIntervalHit = leftRcSAInts.front();
              for (OffsetT sa_i = saIntervalHit.begin; sa_i != saIntervalHit.end; ++sa_i) {
                  auto globalPos = SA[sa_i];
                  auto txpID = sidx->transcriptAtPosition(globalPos);
                  // the offset into this transcript
                  auto pos = globalPos - txpStarts[txpID];
                  int32_t hitPos = pos;                  
                  //int32_t junctPos = 99;//transcripts[txpID].RefLength/2 = 198/2
                  int32_t junctPos = floor(transcripts[txpID].RefLength/2);
//                  std::cout <<std::endl << " txpID " <<txpID << " txName " <<transcripts[txpID].RefName << " txLen " << transcripts[txpID].RefLength << " hitPos "<< std::to_string(hitPos) ; 
                  if ((hitPos < (junctPos-5)) and ((hitPos+saIntervalHit.len) > (junctPos+5))){
                    //FW(0)/RC(1)
                    std::string myJunction=header1+"\t"+ "1" +"\t" + std::to_string(saIntervalHit.queryPos)+"\t"+std::to_string(saIntervalHit.len) +"\t"+std::to_string(hitPos) +"\t"+transcripts[txpID].RefName + "\t"+std::to_string(transcripts[txpID].RefLength);
                    fmt::print(UN_faOutput1.get(), myJunction);
                    fmt::print(UN_faOutput1.get(), "\n");
                    isSplit=true;
                    //std::cout <<std::endl  << myJunction  << std::endl;

                  }
              }
            }


          }



            upperBoundHits += (jointHits.size() > 0);

            // If the read mapped to > maxReadOccs places, discard it
            if (jointHits.size() > sfOpts.maxReadOccs ) { jointHits.clear(); }
            //Circall - unmapped reads then export
            /*
            if (jointHits.size() == 0) {
              fmt::print(UN_faOutput1.get(), header1);
              fmt::print(UN_faOutput1.get(), "\n");
              fmt::print(UN_faOutput1.get(), j->data[i].seq);
              fmt::print(UN_faOutput1.get(), "\n");
            }
            */
            //Circall: just update equivalence classes if there are split reads
            //if (jointHits.size() > 0) {
            if (jointHits.size() > 0 and isSplit) {
                int32_t fwAll = 0;
                int32_t fwCompat = 0;
                int32_t rcAll = 0;
                int32_t rcCompat = 0;

                double auxSumAll = 0.0;
                double auxSumCompat = 0.0;

                bool needBiasSample = sfOpts.biasCorrect;

                for (auto& h : jointHits) {
                    auto transcriptID = h.transcriptID();
                    auto& txp = transcripts[transcriptID];

                    int32_t pos = static_cast<int32_t>(h.pos);
                    auto dir = sailfish::utils::boolToDirection(h.fwd);

                    // Note: sidx is a pointer to type IndexT, not RapMapSAIndex!

                    // If bias correction is turned on, and we haven't sampled a mapping
                    // for this read yet, and we haven't collected the required number of
                    // samples overall.
                    if(needBiasSample and sfOpts.numBiasSamples > 0){
                        // the "start" position is the leftmost position if
                        // we hit the forward strand, and the leftmost
                        // position + the read length if we hit the reverse complement
                        int32_t startPos = h.fwd ? pos : pos + h.readLen;

                        if (startPos > 0 and startPos < txp.RefLength) {
                            /*
                            const char* txpStart = txomeStr + sidx->txpOffsets[h.tid];
                            const char* readStart = txpStart + startPos; // is this correct?
                            const char* txpEnd = txpStart + sidx->txpLens[h.tid]; //??
                            */

                            const char* txpStart = txp.Sequence();
                            const char* readStart = txpStart + startPos; // is this correct?
                            const char* txpEnd = txpStart+ txp.RefLength;
                            bool success = readBias.update(txpStart, readStart, txpEnd, dir);
                            if (success) {
                                sfOpts.numBiasSamples -= 1;
                                needBiasSample = false;
                            }
                        }
                    }

                    // True if the read is compatible with the
                    // expected library type; false otherwise.
                    bool compat = ignoreCompat;
                    if (!compat) {
                        compat = sailfish::utils::compatibleHit(
                                expectedLibType, pos,
                                h.fwd, h.mateStatus);
                    }

                    if (compat) {
                        haveCompat = true;
                        txpIDsCompat.push_back(transcriptID);
                        auxProbsCompat.push_back(1.0);
                        auxSumCompat += 1.0;
                        if (h.fwd) { fwCompat++; } else { rcCompat++; }
                    }
                    if (!haveCompat and !enforceCompat) {
                        txpIDsAll.push_back(transcriptID);
                        auxProbsAll.push_back(1.0);
                        auxSumAll += 1.0;
                        if (h.fwd) { fwAll++; } else { rcAll++; }
                    }
        }

                // If we have compatible hits, only use those
                if (haveCompat) {
                    if (txpIDsCompat.size() > 0) {
                        mappedFrag = true;
                        TranscriptGroup tg(txpIDsCompat);
                        eqBuilder.addGroup(std::move(tg), auxProbsCompat);
                        readExp.addNumFwd(fwCompat);
                        readExp.addNumRC(rcCompat);
                    }
                } else {
                    if (txpIDsAll.size() > 0) {
                        // Otherwise, consider all hits.
                        mappedFrag = true;
                        TranscriptGroup tg(txpIDsAll);
                        eqBuilder.addGroup(std::move(tg), auxProbsAll);
                        readExp.addNumFwd(fwAll);
                        readExp.addNumRC(rcAll);
                    }
                }
            }

            validHits += (mappedFrag) ? 1 : 0;
            totalHits += jointHits.size();
            locRead++;
            ++numObservedFragments;
            if (numObservedFragments % 500000 == 0) {
                iomutex.lock();
                fmt::print(stderr, "\033[A\r\rprocessed {} fragments\n", numObservedFragments);
                fmt::print(stderr, "hits: {}, hits per frag (may not be concordant):  {}",
                        totalHits,
                        totalHits / static_cast<float>(prevObservedFrags));

                iomutex.unlock();
            }

        } // end for i < j->nb_filled

        prevObservedFrags = numObservedFragments;
    }
}

std::vector<double> getNormalFragLengthDist(
        const SailfishOpts& sfOpts) {

    std::vector<double> correctionFactors(sfOpts.maxFragLen, 0.0);
    auto maxLen = sfOpts.maxFragLen;
    auto mean = sfOpts.fragLenDistPriorMean;
    auto sd = sfOpts.fragLenDistPriorSD;

    auto kernel = [mean, sd](double p) -> double {
        double invStd = 1.0 / sd;
        double x = invStd * (p - mean);
        return std::exp(-0.5 * x * x) * invStd;
    };

    double cumulativeMass{0.0};
    double cumulativeDensity{0.0};
    for (size_t i = 0; i < sfOpts.maxFragLen; ++i) {
        auto d = kernel(static_cast<double>(i));
        cumulativeMass += i * d;
        cumulativeDensity += d;
        if (cumulativeDensity > 0) {
            correctionFactors[i] = cumulativeMass / cumulativeDensity;
        }
    }
    return correctionFactors;
}

std::vector<int32_t> getNormalFragLengthCounts(
        const SailfishOpts& sfOpts) {

    std::vector<int> dist(sfOpts.maxFragLen, 0);
    int32_t totalCount = sfOpts.numFragSamples;
    auto maxLen = sfOpts.maxFragLen;
    auto mean = sfOpts.fragLenDistPriorMean;
    auto sd = sfOpts.fragLenDistPriorSD;

    auto kernel = [mean, sd](double p) -> double {
        double invStd = 1.0 / sd;
        double x = invStd * (p - mean);
        return std::exp(-0.5 * x * x) * invStd;
    };

    double totalMass{0.0};
    for (size_t i = 0; i < sfOpts.maxFragLen; ++i) {
        totalMass += kernel(static_cast<double>(i));
    }

    double currentDensity{0.0};
    if (totalMass > 0) {
        for (size_t i = 0; i < sfOpts.maxFragLen; ++i) {
            currentDensity = kernel(static_cast<double>(i));
            dist[i] = static_cast<int>(
                    std::round(currentDensity * totalCount / totalMass));
        }
    }
    return dist;
}


void setEffectiveLengthsDirect(ReadExperiment& readExp,
        const SailfishOpts& sfOpts) {
        auto& transcripts = readExp.transcripts();
        for(size_t txpID = 0; txpID < transcripts.size(); ++txpID) {
            auto& txp = transcripts[txpID];
            double refLen = txp.RefLength;
            txp.EffectiveLength = txp.RefLength;
        }
}

void computeEmpiricalEffectiveLengths(
        const SailfishOpts& sfOpts,
        std::vector<Transcript>& transcripts,
        std::map<uint32_t, uint32_t>& jointMap) {
            std::vector<uint32_t> vals;
            std::vector<uint32_t> multiplicities;

            vals.reserve(jointMap.size());
            multiplicities.reserve(jointMap.size());

            //Circall
            //Extract information of fragment length
            string fragmentDistfn ="fragmentDist.txt";
            boost::filesystem::path fragmentDistfnpath = sfOpts.outputDirectory  / fragmentDistfn;
            std::unique_ptr<std::FILE, int (*)(std::FILE *)> fragmentDistOutput(std::fopen(fragmentDistfnpath.c_str(), "a"), std::fclose); 
            //Circall

            for (auto& kv : jointMap) {
                vals.push_back(kv.first);
                multiplicities.push_back(kv.second);
                //Circall
                std::string mytext=std::to_string(kv.first)+ "\t"+std::to_string(kv.second)+"\n";
                fmt::print(fragmentDistOutput.get(), mytext);
                //Circall
            }

            sfOpts.jointLog->info("Building empirical fragment length distribution");
            EmpiricalDistribution empDist(vals, multiplicities);
            sfOpts.jointLog->info("finished building empirical fragment length distribution");
            using BlockedIndexRange =  tbb::blocked_range<size_t>;

            tbb::task_scheduler_init tbbScheduler(sfOpts.numThreads);

            sfOpts.jointLog->info("Estimating effective lengths");
            sfOpts.jointLog->info("Emp. dist min = {}, Emp. dist max = {}",
                                  empDist.minValue(), empDist.maxValue());

            tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts.size())),
                [&transcripts, &empDist](const BlockedIndexRange& range) -> void {
                    for (auto txpID : boost::irange(range.begin(), range.end())) {
                        auto& txp = transcripts[txpID];
                       /**
                          *  NOTE: Adopted from "est_effective_length" at
                          *  (https://github.com/adarob/eXpress/blob/master/src/targets.cpp)
                          *  originally written by Adam Roberts.
                          */
                        uint32_t minVal = empDist.minValue();
                        uint32_t maxVal = empDist.maxValue();
                        bool validDistSupport = (maxVal > minVal);
                        double refLen = txp.RefLength;
                        double effectiveLength = 0.0;
                        for (size_t l = minVal; l <= std::min(txp.RefLength, maxVal); ++l) {
                            effectiveLength += empDist.pdf(l) * (txp.RefLength - l + 1.0);
                        }
                        if (effectiveLength < 1.0) { 
                            txp.EffectiveLength = txp.RefLength;
                        } else {
                            txp.EffectiveLength = effectiveLength;
                        }
                    }
                });
}

std::vector<double> correctionFactorsFromCounts(
        const SailfishOpts& sfOpts,
        std::map<uint32_t, uint32_t>& jointMap) {
    auto maxLen = sfOpts.maxFragLen;

    std::vector<double> correctionFactors(maxLen, 0.0);
    std::vector<double> vals(maxLen, 0.0);
    std::vector<uint32_t> multiplicities(maxLen, 0);

    auto valIt = jointMap.find(0);
    if (valIt != jointMap.end()) {
        multiplicities[0] = valIt->second;
    } else {
        multiplicities[0] = 0;
    }

    sfOpts.jointLog->info(
            "Computing effective length factors --- max length = {}",
            maxLen);

    uint32_t v{0};
    for (size_t i = 1; i < maxLen; ++i) {
        valIt = jointMap.find(i);
        if (valIt == jointMap.end()) {
            v = 0;
        } else {
            v = valIt->second;
        }
        vals[i] = static_cast<double>(v * i) + vals[i-1];
        multiplicities[i] = v + multiplicities[i-1];
        if (multiplicities[i] > 0) {
            correctionFactors[i] = vals[i] / static_cast<double>(multiplicities[i]);
        }
    }
    sfOpts.jointLog->info("finished computing effective length factors");
    sfOpts.jointLog->info("mean fragment length = {}", correctionFactors[maxLen-1]);

    return correctionFactors;
}

void computeSmoothedEffectiveLengths(
        const SailfishOpts& sfOpts,
        std::vector<Transcript>& transcripts,
        std::vector<double>& correctionFactors) {

            auto maxLen = sfOpts.maxFragLen;
            using BlockedIndexRange =  tbb::blocked_range<size_t>;
            tbb::task_scheduler_init tbbScheduler(sfOpts.numThreads);
            sfOpts.jointLog->info("Estimating effective lengths");

            tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts.size())),
                [&transcripts, &correctionFactors, maxLen](const BlockedIndexRange& range) -> void {

                    for (auto txpID : boost::irange(range.begin(), range.end())) {
                        auto& txp = transcripts[txpID];
                        auto origLen = txp.RefLength;
                        double correctionFactor = (origLen >= maxLen) ?
                                                  correctionFactors[maxLen-1] :
                                                  correctionFactors[origLen];

                        double effLen = static_cast<double>(txp.RefLength) -
                                        correctionFactor + 1.0;
                        if (effLen < 1.0) {
                            effLen = static_cast<double>(origLen);
                        }

                        txp.EffectiveLength = effLen;
                    }
                });
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}


void quasiMapReads(
        ReadExperiment& readExp,
        SailfishOpts& sfOpts,
        std::mutex& iomutex){

    std::vector<std::thread> threads;
    auto& rl = readExp.readLibraries().front();
    rl.checkValid();

    auto numThreads = sfOpts.numThreads;

    std::unique_ptr<paired_parser> pairedParserPtr{nullptr};
    std::unique_ptr<single_parser> singleParserPtr{nullptr};

    // Remember the fragment lengths that we see in each thread
    //std::vector<FragLengthCountMap> flMaps(numThreads);
    FragLengthCountMap flMap(sfOpts.maxFragLen, 0);

    bool largeIndex = readExp.getIndex()->is64BitQuasi();
    bool perfectHashIndex = readExp.getIndex()->isPerfectHashQuasi();

    // If the read library is paired-end
    // ------ Paired-end --------
    if (rl.format().type == ReadType::PAIRED_END) {

        if (rl.mates1().size() != rl.mates2().size()) {
            sfOpts.jointLog->error("The number of provided files for "
                    "-1 and -2 must be the same!");
            sfOpts.jointLog->flush();
            spdlog::drop_all();
            std::this_thread::sleep_for(std::chrono::seconds(1));
            std::exit(1);
        }

        size_t numFiles = rl.mates1().size() + rl.mates2().size();
        char** pairFileList = new char*[numFiles];
        //std::vector<std::ifstream*> pairFileList(numFiles);
        //pairFileList.reserve(numFiles);
        for (size_t i = 0; i < rl.mates1().size(); ++i) {
            pairFileList[2*i] = const_cast<char*>(rl.mates1()[i].c_str());
            pairFileList[2*i+1] = const_cast<char*>(rl.mates2()[i].c_str());
            //pairFileList[2*i] = new std::ifstream(rl.mates1()[i]);
            //pairFileList[2*i+1] = new std::ifstream(rl.mates2()[i]);
        }

        size_t maxReadGroup{readGroupSize}; // Number of reads in each "job"
        size_t concurrentFile{2}; // Number of files to read simultaneously
        pairedParserPtr.reset(new
                paired_parser(4 * numThreads, maxReadGroup,
                    concurrentFile,
                    pairFileList, pairFileList+numFiles));
                    //pairFileList.begin(), pairFileList.end()));

        std::atomic<int32_t> remainingFLOps{sfOpts.numFragSamples};

        for(int i = 0; i < numThreads; ++i)  {
            // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
            // change value before the lambda below is evaluated --- crazy!

            // if we have a 64-bit index
            if (largeIndex) {
        if (perfectHashIndex) {
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int64_t, PerfectHash<int64_t>>>(
                            pairedParserPtr.get(),
                            readExp.getIndex()->quasiIndexPerfectHash64(),
                            readExp,
                            rl,
                            sfOpts,
                            flMap,
                            remainingFLOps,
                            iomutex);
                };
                threads.emplace_back(threadFun);
        } else {
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int64_t, DenseHash<int64_t>>>(
                            pairedParserPtr.get(),
                            readExp.getIndex()->quasiIndex64(),
                            readExp,
                            rl,
                            sfOpts,
                            flMap,
                            remainingFLOps,
                            iomutex);
                };
                threads.emplace_back(threadFun);
        }
            } else { // 32-bit
        
        if (perfectHashIndex) {
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int32_t, PerfectHash<int32_t>>>(
                            pairedParserPtr.get(),
                            readExp.getIndex()->quasiIndexPerfectHash32(),
                            readExp,
                            rl,
                            sfOpts,
                            flMap,
                            remainingFLOps,
                            iomutex);
                };
    threads.emplace_back(threadFun);
        } else {
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int32_t, DenseHash<int32_t>>>(
                            pairedParserPtr.get(),
                            readExp.getIndex()->quasiIndex32(),
                            readExp,
                            rl,
                            sfOpts,
                            flMap,
                            remainingFLOps,
                            iomutex);
                };
    threads.emplace_back(threadFun);
        }

            } // end if (largeIndex) 
        }

        // Join the threads and collect the results from the count maps
        size_t totalObs{0};
        std::map<uint32_t, uint32_t> jointMap;

  // join all the worker threads
        for(int i = 0; i < numThreads; ++i) { threads[i].join(); }

        for (size_t i = 0; i < flMap.size(); ++i) {
            jointMap[i] = flMap[i];
            totalObs += flMap[i];
        }

  // we need an extra newline here.
  fmt::print(stderr, "\n");

        sfOpts.jointLog->info("Gathered fragment lengths from all threads");
        //sfOpts.jointLog->info("total number of mapped fragments > {} is {}", sfOpts.maxFragLen, sfOpts.fragsTooLong.load());


        /** If we have a sufficient number of observations for the empirical
         *  distribution, then use that --- otherwise use the provided prior
         *  mean fragment length.
         **/
        // Note: if "noEffectiveLengthCorrection" is set, so that these values
        // won't matter anyway, then don't bother computing this "expensive"
        // version.
        if (sfOpts.noEffectiveLengthCorrection) {
            setEffectiveLengthsDirect(readExp, sfOpts);
        } else {
            // We didn't have sufficient observations, use the provided
            // values
            if (remainingFLOps > 0) {
                sfOpts.jointLog->warn("Sailfish saw fewer then {} uniquely mapped reads "
                        "so {} will be used as the mean fragment length and {} as "
                        "the standard deviation for effective length correction",
                        sfOpts.numFragSamples,
                        sfOpts.fragLenDistPriorMean,
                        sfOpts.fragLenDistPriorSD);
                // Set the fragment length distribution in the ReadExperiment
                readExp.setFragLengthDist(getNormalFragLengthCounts(sfOpts));
                auto correctionFactors = getNormalFragLengthDist(sfOpts);
                computeSmoothedEffectiveLengths(sfOpts, readExp.transcripts(), correctionFactors);
            } else {
                // Set the fragment length distribution in the ReadExperiment
                std::vector<int32_t> fld(flMap.size(), 0);
                for (size_t i = 0; i < flMap.size(); ++i) {
                    fld[i] = static_cast<int32_t>(flMap[i]);
                }
                readExp.setFragLengthDist(fld);

                if (sfOpts.simplifiedLengthCorrection) {
                    auto correctionFactors = correctionFactorsFromCounts(sfOpts, jointMap);
                    computeSmoothedEffectiveLengths(sfOpts, readExp.transcripts(), correctionFactors);
                } else {
                    computeEmpiricalEffectiveLengths(sfOpts, readExp.transcripts(), jointMap);
                }
            }
        }
    } // ------ Single-end --------
    else if (rl.format().type == ReadType::SINGLE_END) {

        char* readFiles[] = { const_cast<char*>(rl.unmated().front().c_str()) };
        size_t maxReadGroup{readGroupSize}; // Number of files to read simultaneously
        size_t concurrentFile{1}; // Number of reads in each "job"
        stream_manager streams( rl.unmated().begin(),
                rl.unmated().end(), concurrentFile);

        singleParserPtr.reset(new single_parser(4 * numThreads,
                    maxReadGroup,
                    concurrentFile,
                    streams));

        for(int i = 0; i < numThreads; ++i)  {
            // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
            // change value before the lambda below is evaluated --- crazy!
            if (largeIndex) {
        
        if (perfectHashIndex) {
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int64_t, PerfectHash<int64_t>>>(
                            singleParserPtr.get(),
                            readExp.getIndex()->quasiIndexPerfectHash64(),
                            readExp,
                            rl,
                            sfOpts,
                            iomutex);
                };
                threads.emplace_back(threadFun);
        } else { // dense hash index
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int64_t, DenseHash<int64_t>>>(
                            singleParserPtr.get(),
                            readExp.getIndex()->quasiIndex64(),
                            readExp,
                            rl,
                            sfOpts,
                            iomutex);
                };
                threads.emplace_back(threadFun);
        }

            } else { // 32-bit
        
        if (perfectHashIndex) {
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int32_t, PerfectHash<int32_t>>>(
                            singleParserPtr.get(),
                            readExp.getIndex()->quasiIndexPerfectHash32(),
                            readExp,
                            rl,
                            sfOpts,
                            iomutex);
                };
    threads.emplace_back(threadFun);
        } else { // dense hash index
                auto threadFun = [&,i]() -> void {
      processReadsQuasi<RapMapSAIndex<int32_t, DenseHash<int32_t>>>(
                            singleParserPtr.get(),
                            readExp.getIndex()->quasiIndex32(),
                            readExp,
                            rl,
                            sfOpts,
                            iomutex);
                };
    threads.emplace_back(threadFun);
        }

            } // end if (largeIndex)
        }
        for(int i = 0; i < numThreads; ++i) { threads[i].join(); }
        if (sfOpts.noEffectiveLengthCorrection) {
            setEffectiveLengthsDirect(readExp, sfOpts);
        } else {
            // Set the fragment length distribution in the ReadExperiment
            readExp.setFragLengthDist(getNormalFragLengthCounts(sfOpts));

            auto correctionFactors = getNormalFragLengthDist(sfOpts);
            computeSmoothedEffectiveLengths(sfOpts, readExp.transcripts(), correctionFactors);
        }
    } // ------ END Single-end --------
}



int main(int argc, char* argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool biasCorrect{false};
    SailfishOpts sopt;
    sopt.numThreads = std::thread::hardware_concurrency();
    sopt.allowOrphans = true;
    int32_t numBiasSamples{0};

    vector<string> unmatedReadFiles;
    vector<string> mate1ReadFiles;
    vector<string> mate2ReadFiles;
    string txpAggregationKey;

    bool discardOrphans = false;
    po::options_description generic("\n"
            "basic options");
    generic.add_options()
        ("version,v", "print version string")
        ("help,h", "produce help message")
        ("index,i", po::value<string>()->required(), "Sailfish index")
//        ("libType,l", po::value<std::string>()->required(), "Format string describing the library type")
        ("libType,l", po::value<std::string>()->default_value("IU"), "Format string describing the library type")
        ("unmatedReads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
         "List of files containing unmated reads of (e.g. single-end reads)")
        ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
         "File containing the #1 mates")
        ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
         "File containing the #2 mates")
        ("threads,p", po::value<uint32_t>(&(sopt.numThreads))->default_value(sopt.numThreads), "The number of threads to use concurrently.")
        ("output,o", po::value<std::string>()->required(), "Output quantification file.")
        ("geneMap,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided "
         "Sailfish will output both quant.sf and quant.genes.sf files, where the latter "
         "contains aggregated gene-level abundance estimates.  The transcript to gene mapping "
         "should be provided as either a GTF file, or a in a simple tab-delimited format "
         "where each line contains the name of a transcript and the gene to which it belongs "
         "separated by a tab.  The extension of the file is used to determine how the file "
         "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF "
         "format; files with any other extension are assumed to be in the simple format.")
       ("biasCorrect", po::value(&(sopt.biasCorrect))->zero_tokens(), "Perform sequence-specific bias correction")
       ("gcBiasCorrect", po::value(&(sopt.gcBiasCorrect))->zero_tokens(), "[experimental] Perform fragment GC bias correction");



    po::options_description advanced("\n"
            "advanced options");
    advanced.add_options()
        ("auxDir", po::value<std::string>(&(sopt.auxDir))->default_value("aux"), "The sub-directory of the quantification directory where auxiliary information "
          "e.g. bootstraps, bias parameters, etc. will be written.")
        ("dumpEq", po::bool_switch(&(sopt.dumpEq))->default_value(false), "Dump the equivalence class counts "
            "that were computed during quasi-mapping")
        ("gcSizeSamp", po::value<std::uint32_t>(&(sopt.gcSampFactor))->default_value(1), "The value by which to down-sample transcripts when representing the "
             "GC content.  Larger values will reduce memory usage, but may decrease the fidelity of bias modeling results.")
        ("gcSpeedSamp", po::value<std::uint32_t>(&(sopt.pdfSampFactor))->default_value(1), "The value at which the fragment length PMF is down-sampled "
             "when evaluating GC fragment bias.  Larger values speed up effective length correction, but may decrease the fidelity of bias modeling results.")
        ("strictIntersect", po::bool_switch(&(sopt.strictIntersect))->default_value(false), "Modifies how orphans are "
            "assigned.  When this flag is set, if the intersection of the quasi-mappings for the left and right "
            "is empty, then all mappings for the left and all mappings for the right read are reported as orphaned "
            "quasi-mappings")
        ("simplifiedLengthCorrection", po::bool_switch(&(sopt.simplifiedLengthCorrection))->default_value(false), "Use a \"simplfied\" "
            "effective length correction approach, rather than convolving the FLD with the "
            "characteristic function over each transcript.")
        ("maxFragLen", po::value<uint32_t>(&(sopt.maxFragLen))->default_value(1000), "The maximum length of a fragment to consider when "
            "building the empirical fragment length distribution")
        //("readEqClasses", po::value<std::string>(&eqClassFile), "Read equivalence classes in directly")
        ("txpAggregationKey", po::value<std::string>(&txpAggregationKey)->default_value("gene_id"), "When generating the gene-level estimates, "
            "use the provided key for aggregating transcripts.  The default is the \"gene_id\" field, but other fields (e.g. \"gene_name\") might "
            "be useful depending on the specifics of the annotation being used.  Note: this option only affects aggregation when using a "
            "GTF annotation; not an annotation in \"simple\" format.")
        ("ignoreLibCompat", po::bool_switch(&(sopt.ignoreLibCompat))->default_value(false), "Disables "
             "strand-aware processing completely.  All hits are considered \"valid\".")
        ("enforceLibCompat", po::bool_switch(&(sopt.enforceLibCompat))->default_value(false), "Enforces "
             "\"strict\" library compatibility.  Fragments that map in a manner other than what is "
             "specified by the expected library type will be discarded, even if there are no mappings that "
             "agree with the expected library type.")
        ("allowDovetail", po::bool_switch(&(sopt.allowDovetail))->default_value(false), "Allow "
             "paired-end reads from the same fragment to \"dovetail\", such that the ends "
             "of the mapped reads can extend past each other.")
        ("discardOrphans", po::bool_switch(&discardOrphans)->default_value(false), "This option will discard orphaned fragments.  This only "
            "has an effect on paired-end input, but enabling this option will discard, rather than count, any reads where only one of the paired "
            "fragments maps to a transcript.")
        ("noBiasLengthThreshold", po::bool_switch(&(sopt.noBiasLengthThreshold))->default_value(false), "[experimental] : "
                        "If this option is enabled, then bias correction will be allowed to estimate effective lengths "
                        "shorter than the approximate mean fragment length")
        ("numBiasSamples", po::value<int32_t>(&numBiasSamples)->default_value(1000000),
            "Number of fragment mappings to use when learning the sequence-specific bias model.")
        ("numFragSamples", po::value<int32_t>(&(sopt.numFragSamples))->default_value(10000),
            "Number of fragments from unique alignments to sample when building the fragment "
            "length distribution")
        ("fldMean", po::value<size_t>(&(sopt.fragLenDistPriorMean))->default_value(200),
            "If single end reads are being used for quantification, or there are an insufficient "
            "number of uniquely mapping reads when performing paired-end quantification to estimate "
            "the empirical fragment length distribution, then use this value to calculate effective lengths.")
        ("fldSD" , po::value<size_t>(&(sopt.fragLenDistPriorSD))->default_value(80),
            "The standard deviation used in the fragment length distribution for single-end quantification or "
            "when an empirical distribution cannot be learned.")
        ("maxReadOcc,w", po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(200), "Reads \"mapping\" to more than this many places won't be considered.")
        ("noEffectiveLengthCorrection", po::bool_switch(&(sopt.noEffectiveLengthCorrection))->default_value(false), "Disables "
         "effective length correction when computing the probability that a fragment was generated "
         "from a transcript.  If this flag is passed in, the fragment length distribution is not taken "
         "into account when computing this probability.")
        ("useVBOpt", po::bool_switch(&(sopt.useVBOpt))->default_value(false), "Use the Variational Bayesian EM rather than the "
          "traditional EM algorithm to estimate transcript abundances.")
        ("numGibbsSamples", po::value<uint32_t>(&(sopt.numGibbsSamples))->default_value(0), "[*super*-experimental]: Number of Gibbs sampling rounds to "
            "perform.")
        ("numBootstraps", po::value<uint32_t>(&(sopt.numBootstraps))->default_value(0), "[*super*-experimental]: Number of bootstrap samples to generate. Note: "
            "This is mutually exclusive with Gibbs sampling.");

    po::options_description all("sailfish quant options");
    all.add(generic).add(advanced);

    po::options_description visible("sailfish quant options");
    visible.add(generic).add(advanced);

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(all).run();

        po::store(orderedOptions, vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
                Circall_bsj
                ==========
                To perform quasi-mapping to back-splicing-junction database
                This tool was adapted from sailfish
                )";
            std::cout << hstring << std::endl;
            std::cout << visible << std::endl;
            std::exit(1);
        }

        po::notify(vm);

        if (discardOrphans) {
            sopt.allowOrphans = false;
        }
/*
        std::stringstream commentStream;
        commentStream << "# sailfish (quasi-mapping-based) v" << sailfish::version << "\n";
        commentStream << "# [ program ] => sailfish \n";
        commentStream << "# [ command ] => quant \n";
        for (auto& opt : orderedOptions.options) {
            commentStream << "# [ " << opt.string_key << " ] => {";
            for (auto& val : opt.value) {
                commentStream << " " << val;
            }
            commentStream << " }\n";
        }
        std::string commentString = commentStream.str();
        fmt::print(stderr, "{}", commentString);
*/
        // Set the atomic variable numBiasSamples from the local version.
        sopt.numBiasSamples.store(numBiasSamples);

        // Get the time at the start of the run
        std::time_t result = std::time(NULL);
        std::string runStartTime(std::asctime(std::localtime(&result)));
        runStartTime.pop_back(); // remove the newline

        // Verify the geneMap before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("geneMap")) {
            // Make sure the provided file exists
            geneMapPath = vm["geneMap"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                std::cerr << "Could not find transcript <=> gene map file " << geneMapPath << "\n";
                std::cerr << "Exiting now: please either omit the \'geneMap\' option or provide a valid file\n";
                return 1;
            }
        }

        bfs::path outputDirectory(vm["output"].as<std::string>());
        bfs::create_directories(outputDirectory);
        if (!(bfs::exists(outputDirectory) and bfs::is_directory(outputDirectory))) {
            std::cerr << "Couldn't create output directory " << outputDirectory << "\n";
            std::cerr << "exiting\n";
            return 1;
        }

        bfs::path indexDirectory(vm["index"].as<string>());
        bfs::path logDirectory = outputDirectory / "LogDir";

        sopt.indexDirectory = indexDirectory;
        sopt.outputDirectory = outputDirectory;

        // Create the logger and the logging directory
        bfs::create_directories(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "Logs will be written to " << logDirectory.string() << "\n";

        bfs::path logPath = logDirectory / "Circall.log";
        // must be a power-of-two
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        std::ofstream logFile(logPath.string());
        if (!logFile.good()) {
            std::cerr << "[WARNING]: Could not open log file --- this seems suspicious!\n";
        }

        auto fileSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(logFile);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("stderrLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        sopt.jointLog = jointLog;
        sopt.fileLog = fileLog;

/*
        // Write out information about the command / run
        {
            bfs::path cmdInfoPath = outputDirectory / "cmd_info.json";
            std::ofstream os(cmdInfoPath.string());
            {
                cereal::JSONOutputArchive oa(os);
                oa(cereal::make_nvp("sf_version", std::string(sailfish::version)));
                for (auto& opt : orderedOptions.options) {
                    if (opt.value.size() == 1) {
                        oa(cereal::make_nvp(opt.string_key, opt.value.front()));
                    } else {
                        oa(cereal::make_nvp(opt.string_key, opt.value));
                    }
                }
            }
            os.close();
        }
*/
        jointLog->info("parsing read library format");

        if (sopt.numGibbsSamples > 0 and sopt.numBootstraps > 0) {
            jointLog->error("You cannot perform both Gibbs sampling and bootstrapping. "
                            "Please choose one.");
            jointLog->flush();
            spdlog::drop_all();
            return 1;
        }

        vector<ReadLibrary> readLibraries = sailfish::utils::extractReadLibraries(orderedOptions);

        // Verify that no inconsistent options were provided
        {
    
          if (sopt.biasCorrect and sopt.gcBiasCorrect) {
            sopt.jointLog->warn("\n\n"
                 "===================================\n"
                 "Enabling both sequence-specific and fragment GC bias correction "
                 "simultaneously is still experimental, we currently recommend you enable only one for a given sample.\n"
                 "===================================\n\n");
          }
    
          if (sopt.gcBiasCorrect) {
            for (auto& rl : readLibraries) {
              // We can't use fragment GC correction with single
              // end reads yet.
              if (rl.format().type == ReadType::SINGLE_END) {
                jointLog->warn("Fragment GC bias correction is currently "
                    "only implemented for paired-end libraries. "
                    "It is being disabled");
                sopt.gcBiasCorrect = false;
                break;
              }
            }
          }
        } // Done verifying options

        SailfishIndexVersionInfo versionInfo;
        boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
        versionInfo.load(versionPath);

        ReadExperiment experiment(readLibraries, indexDirectory, sopt);
        // end parameter validation

        // This will be the class in charge of maintaining our
        // rich equivalence classes
        experiment.equivalenceClassBuilder().start();
        std::mutex ioMutex;
        fmt::print(stderr, "\n\n");
        quasiMapReads(experiment, sopt, ioMutex);
        fmt::print(stderr, "Done Quasi-Mapping \n\n");
        experiment.equivalenceClassBuilder().finish();

        // Now that we have our reads mapped and our equivalence
        // classes, iterate the abundance estimates to convergence.
        ExportFeq feqExporter;
        jointLog->info("Start to export equivalence classes:\n");
        //for Circall
        boost::filesystem::path feqOutputPath = outputDirectory;
        bool expSuccess = feqExporter.writeFeq(experiment, sopt,feqOutputPath);//for Circall
        if (!expSuccess) {
            jointLog->error("There are errors when exporting equivalence classes.\n"
                            "Please contact us in our website.\n");
            return 1;
        }
        jointLog->info("Finished exporting equivalence classes");


        jointLog->flush();
        spdlog::drop_all();
        logFile.close();

    } catch (po::error &e) {
        std::cerr << "Exception: [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception: [" << e.what() << "]\n";
        std::cerr << argv[0] << " tool was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " Circall --help\nExiting.\n";
        std::exit(1);
    }
    return 0;
}


