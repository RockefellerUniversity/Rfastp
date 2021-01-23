#include <Rcpp.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "GENCORE/cmdline.h"
#include "GENCORE/common.h"
#include <sstream>
#include "GENCORE/util.h"
#include "GENCORE/gencore.h"
#include "GENCORE/options.h"
#include "GENCORE/reference.h"
#include "GENCORE/unittest.h"

using namespace std;
//using namespace Rcpp;

//string gencore_command;

// [[Rcpp::export]]
int runGencore(std::string inBam="",
               std::string outBam="",
               std::string refFile="",
               std::string bedFile="",
               std::string umiPrefix="",
               int numSupportingReads=1,
               int majorBaseScore=6,
               double majorBaseRatio=0.8,
               int quitAfterContig=0,
               int highQual=30,
               int moderateQual=20,
               int lowQual=15,
               int coverageSampling=10000,
               bool debug=false,
           bool verbose=true
               ){


//
  GENCORE::Options opt;
  opt.input = inBam;
  opt.output = outBam;
  opt.refFile = refFile;
  opt.bedFile = bedFile;
  opt.umiPrefix = umiPrefix;
  opt.clusterSizeReq = numSupportingReads;
  opt.baseScoreReq = majorBaseScore;
  opt.scorePercentReq = majorBaseRatio;
  opt.maxContig = quitAfterContig;
  opt.highQuality = highQual;
  opt.moderateQuality = moderateQual;
  opt.lowQuality = lowQual;
  opt.coverageStep = coverageSampling;
  opt.debug = debug;
  opt.jsonFile = outBam + ".json";
  opt.htmlFile = outBam + ".html";
  opt.verbose = verbose;
//
//   opt.validate();
//
//   time_t t1 = time(NULL);

  // loading reference
  GENCORE::Reference* reference = NULL;
  if(!opt.refFile.empty() && verbose) {
    Rcpp::Rcerr << "loading reference data:" << endl;
    reference = GENCORE::Reference::instance(&opt);
  }


  GENCORE::Gencore gencore(&opt);
  gencore.consensus();
  //
  // // if(reference) {
  // //   delete reference;
  // //   reference=NULL;
  // // }
  //
  // time_t t2 = time(NULL);
  // // cerr << endl << command << endl;
  // cerr << "gencore v" << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds" << endl;

  return 0;
}
