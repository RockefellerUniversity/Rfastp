#include "options.h"
#include "util.h"
#include <Rcpp.h>

namespace GENCORE {
Options::Options(){
    input = "";
    output = "";
    refFile = "";
    bedFile = "";
    umiPrefix = "";
    maxContig = 0;
    bamHeader = NULL;
    properReadsUmiDiffThreshold = 2;
    unproperReadsUmiDiffThreshold = 0;
    debug = false;
    hasBedFile = false;

    scorePercentReq = 0.8;
    clusterSizeReq = 1;

    highQuality = 30;
    moderateQuality = 20;
    lowQuality = 15;

    scoreOfNotOverlapped = 6;
    scoreOfHighQualityMatch = 8;
    scoreOfLowQualityMatch = 7;
    scoreOfBothHighQualityMismatch = 4;
    scoreOfBothModerateQualityMismatch = 3;
    scoreOfBothLowQualityMismatch = 2;
    scoreOfUnbalancedMismatchHighQuality = 5;
    scoreOfUnbalancedMismatchLowQuality = 1;

    baseScoreReq = scoreOfNotOverlapped;
    skipLowComplexityClusterThreshold = 1000;

    reportTitle = "gencore report";

    bedCoverageStep = 10;
    coverageStep = 10000;
}

bool Options::validate() {
    if(input.empty()) {
        Rcpp::stop("input should be specified by --in1");
    } else {
        check_file_valid(input);
    }

    if(ends_with(refFile, ".gz") || ends_with(refFile, ".gz")) {
        Rcpp::Rcerr << "reference fasta file should not be compressed.\nplease unzip "<<refFile<<" and try again."<<endl;
	Rcpp::stop("\n");
        //exit(-1);
    }

    if(scorePercentReq > 1.0) {
        Rcpp::stop("ratio_threshold cannot be greater than 1.0");
    } else if(scorePercentReq < 0.5) {
        Rcpp::stop("ratio_threshold cannot be less than 0.5");
    }

    if(clusterSizeReq > 10) {
        Rcpp::stop("supporting_reads cannot be greater than 10");
    } else if(clusterSizeReq < 1) {
        Rcpp::stop("supporting_reads cannot be less than 1");
    }

    if(baseScoreReq > 10) {
        Rcpp::stop("score_threshold cannot be greater than 10");
    } else if(baseScoreReq < 1) {
        Rcpp::stop("score_threshold cannot be less than 1");
    }

    if(highQuality > 40) {
        Rcpp::stop("high_qual cannot be greater than 40");
    } else if(highQuality < 20) {
        Rcpp::stop("high_qual cannot be less than 20");
    }

    if(moderateQuality > 35) {
        Rcpp::stop("moderate_qual cannot be greater than 35");
    } else if(moderateQuality < 15) {
        Rcpp::stop("moderate_qual cannot be less than 15");
    }

    if(lowQuality > 30) {
        Rcpp::stop("low_qual cannot be greater than 30");
    } else if(lowQuality < 8) {
        Rcpp::stop("low_qual cannot be less than 8");
    }

    if(lowQuality > moderateQuality) {
        Rcpp::stop("low_qual cannot be greater than moderate_qual");
    }

    if(moderateQuality > highQuality) {
        Rcpp::stop("moderate_qual cannot be greater than high_qual");
    }

    return true;
}
}
