#include <Rcpp.h>

#include <stdio.h>
#include "FASTP/fastqreader.h"
#include "FASTP/unittest.h"
#include <time.h>
#include "FASTP/cmdline.h"
#include <sstream>
#include "FASTP/util.h"
#include "FASTP/options.h"
#include "FASTP/processor.h"
#include "FASTP/evaluator.h"

string command = "testcommand";
using namespace std;
//using namespace Rcpp;
mutex logmtx;

//             int compressLevel=4,
//             int readsToProcess=0,
//             bool dontOverwrite=false,
//             int splitPrefixPaddingNum=4,
//             std::string reportTitle="Rfastp Report",


// [[Rcpp::export]]
int runFastp(std::string read1="",
             std::string read2="",
             std::string outputFastq="",
             std::string unpaired="",
             std::string failedOut="",
             bool merge=false,
             std::string mergeOut="",
             bool phred64=false,
             bool interleaved=false,
             bool fixMGIid=false,
             bool adapterTrimming=true,
             std::string adapterSequenceRead1="auto",
             std::string adapterSequenceRead2="auto",
             std::string adapterFasta="",
             int trimFrontRead1=0,
             int trimTailRead1=0,
             int trimFrontRead2=0,
             int trimTailRead2=0,
             int maxLengthRead1=0,
             int maxLengthRead2=0,
             bool forceTrimPolyG=false,
             bool disableTrimPolyG=false,
             int minLengthPolyG=10,
             bool trimPolyX=false,
             int minLengthPolyX=10,
             bool cutLowQualFront=false,
             bool cutLowQualTail=false,
             bool cutSlideWindowRight=false,
             int cutWindowSize=4,
             int cutMeanQual=20,
             int cutFrontWindowSize=4,
             int cutFrontMeanQual=20,
             int cutTailWindowSize=4,
             int cutTailMeanQual=20,
             int cutSlideWindowSize=4,
             int cutSlideWindowQual=20,
             bool qualityFiltering=true,
             int qualityFilterPhred=15,
             int qualityFilterPercent=40,
             int maxNfilter=5,
             int averageQualFilter=0,
             bool lengthFiltering=true,
             int minReadLength=15,
             int maxReadLength=0,
             bool lowComplexityFiltering=false,
             int minComplexity=30,
             std::string index1Filter="",
             std::string index2Filter="",
             int maxIndexMismatch=0,
             bool correctionOverlap=false,
             int minOverlapLength=30,
             int maxOverlapMismatch=5,
             int maxOverlapMismatchPercentage=20,
             bool umi=false,
             bool umiIgnoreSeqNameSpace=false,
             std::string umiLoc="",
             int umiLength=0,
             std::string umiPrefix="",
             int umiSkipBaseLength=0,
             bool overrepresentationAnalysis=false,
             int overrepresentationSampling=20,
             int splitOutput=0,
             int splitByLines=0,
             int thread=2,
             bool verbose=true){
    Options opt;

    // I/O
    opt.in1 = read1;
    opt.in2 = read2;
    opt.out1 = outputFastq + "_R1.fastq.gz";
    opt.out2 = "";
    if (!opt.in2.empty() || interleaved)
        opt.out2 = outputFastq + "_R2.fastq.gz";
    opt.unpaired1 = unpaired;
    opt.unpaired2 = unpaired;
    opt.failedOut = failedOut;
    opt.overlappedOut = mergeOut;
//    int compressLevel=4;
//    int readsToProcess=0;
    opt.compression = 4;
    opt.readsToProcess = 0;
    opt.phred64 = phred64;
    opt.dontOverwrite = false;
    opt.inputFromSTDIN = false;
    opt.outputToSTDOUT = false;
    opt.interleavedInput = interleaved;
    opt.verbose = verbose;
    opt.fixMGI = fixMGIid;

    // merge PE
    opt.merge.enabled = merge;
    opt.merge.out = mergeOut;
    opt.merge.includeUnmerged = false;

    // adapter cutting
    opt.adapter.enabled = adapterTrimming;
    opt.adapter.detectAdapterForPE = true;
    opt.adapter.sequence = adapterSequenceRead1;
    opt.adapter.sequenceR2 = adapterSequenceRead2;
    opt.adapter.fastaFile = adapterFasta;
    if(opt.adapter.sequenceR2=="auto" && !opt.adapter.detectAdapterForPE && opt.adapter.sequence != "auto") {
        opt.adapter.sequenceR2 = opt.adapter.sequence;
    }
    if(!opt.adapter.fastaFile.empty()) {
        opt.loadFastaAdapters();
    }

    // trimming
    opt.trim.front1 = trimFrontRead1;
    opt.trim.tail1 = trimTailRead1;
    opt.trim.maxLen1 = maxLengthRead1;
    opt.trim.front2 = trimFrontRead2;
    opt.trim.tail2 = trimTailRead2;
    opt.trim.maxLen2 = maxLengthRead2;

    // polyG tail trimming
    if( forceTrimPolyG && disableTrimPolyG) {
        Rcpp::stop("You cannot enabled both trim_poly_g and disable_trim_poly_g");
    } else if(forceTrimPolyG) {
        opt.polyGTrim.enabled = true;
    } else if(disableTrimPolyG) {
        opt.polyGTrim.enabled = false;
    }
    opt.polyGTrim.minLen = minLengthPolyG;

    // polyX tail trimming
    if(trimPolyX) {
        opt.polyXTrim.enabled = true;
    }
    opt.polyXTrim.minLen = minLengthPolyX;


    // sliding window cutting by quality
    opt.qualityCut.enabledFront = cutLowQualFront;
    opt.qualityCut.enabledTail = cutLowQualTail;
    opt.qualityCut.enabledRight = cutSlideWindowRight;
    opt.qualityCut.windowSizeShared = cutWindowSize;
    opt.qualityCut.qualityShared = cutMeanQual;

    opt.qualityCut.windowSizeFront = cutFrontWindowSize;
    opt.qualityCut.qualityFront = cutFrontMeanQual;

    opt.qualityCut.windowSizeTail = cutTailWindowSize;
    opt.qualityCut.qualityTail = cutTailMeanQual;

    opt.qualityCut.windowSizeRight = cutSlideWindowSize;
    opt.qualityCut.qualityRight = cutSlideWindowQual;

    // ============= default parameters issues =============
    // raise a warning if cutting option is not enabled but -W/-M is enabled
    //if(!opt.qualityCut.enabledFront && !opt.qualityCut.enabledTail && !opt.qualityCut.enabledRight) {
    //    if(cutMeanQual > 0 )
    //        Rcpp::Rcerr << "WARNING: you specified the options for cutting by quality, but forogt to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    //}

    // quality filtering
    opt.qualfilter.enabled = qualityFiltering;
    opt.qualfilter.qualifiedQual = num2qual(qualityFilterPhred);
    opt.qualfilter.unqualifiedPercentLimit = qualityFilterPercent;
    opt.qualfilter.avgQualReq = averageQualFilter;
    opt.qualfilter.nBaseLimit = maxNfilter;

    // length filtering
    opt.lengthFilter.enabled = lengthFiltering;
    opt.lengthFilter.requiredLength = minReadLength;
    opt.lengthFilter.maxLength = maxReadLength;

    // low complexity filter
    opt.complexityFilter.enabled = lowComplexityFiltering;
    opt.complexityFilter.threshold = (min(100, max(0, minComplexity))) / 100.0;

    // overlap correction
    opt.correction.enabled = correctionOverlap;
    opt.overlapRequire = minOverlapLength;
    opt.overlapDiffLimit = maxOverlapMismatch;
    opt.overlapDiffPercentLimit = maxOverlapMismatchPercentage;

    // threading
    opt.thread = thread;

    // reporting
    opt.jsonFile = outputFastq + ".json";
    opt.htmlFile = outputFastq + ".html";
//    opt.reportTitle = reportTitle;
    std::string reportTitle="Rfastp Report";
    opt.reportTitle = reportTitle;

    // splitting
    opt.split.enabled = splitOutput > 0 || splitByLines > 0;
    int splitPrefixPaddingNum=4;
    opt.split.digits = splitPrefixPaddingNum;
    if( splitOutput > 0  && splitByLines > 0 ) {
        Rcpp::stop("You cannot set both splitting by file number and splitting by file lines, please choose either.");
    }
    if(splitOutput > 0) {
        opt.split.number = splitOutput;
        opt.split.needEvaluation = true;
        opt.split.byFileNumber = true;
    }
    if(splitByLines > 0) {
        long lines = splitByLines;
        if(lines % 4 != 0) {
            Rcpp::stop("Line number should be a multiple of 4");
        }
        opt.split.size = lines / 4; // 4 lines per record
        opt.split.needEvaluation = false;
        opt.split.byFileLines = true;
    }

    // umi
    opt.umi.enabled = umi;
    opt.umi.ignore = umiIgnoreSeqNameSpace;
    opt.umi.length = umiLength;
    opt.umi.prefix = umiPrefix;
    opt.umi.skip = umiSkipBaseLength;
    if(opt.umi.enabled) {
        //string umiLoc = umiLoc;
        str2lower(umiLoc);
        if(umiLoc.empty())
            Rcpp::stop("You've enabled UMI (umi=TRUE), you should specify the UMI location by umiLoc");
        if(umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" && umiLoc != "per_index" && umiLoc != "per_read") {
            Rcpp::stop("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }
        if(!opt.isPaired() && (umiLoc == "index2" || umiLoc == "read2"))
            Rcpp::stop("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if(opt.umi.length == 0 && (umiLoc == "read1" || umiLoc == "read2" ||  umiLoc == "per_read"))
            Rcpp::stop("You specified the UMI location as " + umiLoc + ", but the length is not specified (--umi_len).");
        if(umiLoc == "index1") {
            opt.umi.location = UMI_LOC_INDEX1;
        } else if(umiLoc == "index2") {
            opt.umi.location = UMI_LOC_INDEX2;
        } else if(umiLoc == "read1") {
            opt.umi.location = UMI_LOC_READ1;
        } else if(umiLoc == "read2") {
            opt.umi.location = UMI_LOC_READ2;
        } else if(umiLoc == "per_index") {
            opt.umi.location = UMI_LOC_PER_INDEX;
        } else if(umiLoc == "per_read") {
            opt.umi.location = UMI_LOC_PER_READ;
        }
    }

    // overrepresented sequence analysis
    opt.overRepAnalysis.enabled = overrepresentationAnalysis;
    opt.overRepAnalysis.sampling = overrepresentationSampling;

    // filtering by index
    string blacklist1 = index1Filter;
    string blacklist2 = index2Filter;
    int indexFilterThreshold = maxIndexMismatch;
    opt.initIndexFiltering(blacklist1, blacklist2, indexFilterThreshold);
/*
    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    bool supportEvaluation = !opt.inputFromSTDIN && opt.in1!="/dev/stdin";

    Evaluator eva(&opt);
    if(supportEvaluation) {
        eva.evaluateSeqLen();

        if(opt.overRepAnalysis.enabled)
            eva.evaluateOverRepSeqs();
    }
    long readNum = 0;
*/

    // using evaluator to guess how many reads in total
/*    if(opt.shallDetectAdapter(false)) {
        if(!supportEvaluation)
            cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        else {
            cerr << "Detecting adapter sequence for read1..." << endl;
            string adapt = eva.evalAdapterAndReadNum(readNum, false);
            if(adapt.length() > 60 )
                adapt.resize(0, 60);
            if(adapt.length() > 0 ) {
                opt.adapter.sequence = adapt;
                opt.adapter.detectedAdapter1 = adapt;
            } else {
                cerr << "No adapter detected for read1" << endl;
                opt.adapter.sequence = "";
            }
            cerr << endl;
        }
    }
    if(opt.shallDetectAdapter(true)) {
        if(!supportEvaluation)
            cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        else {
            cerr << "Detecting adapter sequence for read2..." << endl;
            string adapt = eva.evalAdapterAndReadNum(readNum, true);
            if(adapt.length() > 60 )
                adapt.resize(0, 60);
            if(adapt.length() > 0 ) {
                opt.adapter.sequenceR2 = adapt;
                opt.adapter.detectedAdapter2 = adapt;
            } else {
                cerr << "No adapter detected for read2" << endl;
                opt.adapter.sequenceR2 = "";
            }
            cerr << endl;
        }
    }

    opt.validate();
*/
    // using evaluator to guess how many reads in total
/*    if(opt.split.needEvaluation && supportEvaluation) {
        // if readNum is not 0, means it is already evaluated by other functions
        if(readNum == 0) {
            eva.evaluateReadNum(readNum);
        }
        opt.split.size = readNum / opt.split.number;
        // one record per file at least
        if(opt.split.size <= 0) {
            opt.split.size = 1;
            cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
        }
    }
*/
    // using evaluator to check if it's two color system
/*    if(!cmd.exist("trim_poly_g") && !cmd.exist("disable_trim_poly_g") && supportEvaluation) {
        bool twoColorSystem = eva.isTwoColorSystem();
        if(twoColorSystem){
            opt.polyGTrim.enabled = true;
        }
    }
*/
    Processor p(&opt);
    p.process();
/*
    time_t t2 = time(NULL);

    cerr << endl << "JSON report: " << opt.jsonFile << endl;
    cerr << "HTML report: " << opt.htmlFile << endl;
    cerr << endl << command << endl;
    cerr << "fastp v" << FASTP_VER << ", time used: " << (t2)-t1 << " seconds" << endl;
*/
    return 0;
}
