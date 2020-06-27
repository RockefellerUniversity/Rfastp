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
int runGencore(std::string input="-",
               std::string output="-",
               std::string refFile="",
               std::string bedFile="",
               std::string umiPrefix="",
               std::string jsonFile="gencore.json",
               std::string htmlFile="gencore.html",
               int clusterSizeReq=1,
               int baseScoreReq=6,
               double scorePercentReq=0.8,
               int maxContig=0,
               int highQuality=30,
               int moderateQuality=20,
               int lowQuality=15,
               int coverageStep=10000,
               bool debug=false
               ){


//
  GENCORE::Options opt;
  opt.input = input;
  opt.output = output;
  opt.refFile = refFile;
  opt.bedFile = bedFile;
  opt.umiPrefix = umiPrefix;
  opt.clusterSizeReq = clusterSizeReq;
  opt.baseScoreReq = baseScoreReq;
  opt.scorePercentReq = scorePercentReq;
  opt.maxContig = maxContig;
  opt.highQuality = highQuality;
  opt.moderateQuality = moderateQuality;
  opt.lowQuality = lowQuality;
  opt.coverageStep = coverageStep;
  opt.debug = debug;
  opt.jsonFile = jsonFile;
  opt.htmlFile = htmlFile;
//
//   opt.validate();
//
//   time_t t1 = time(NULL);

  // loading reference
  GENCORE::Reference* reference = NULL;
  if(!opt.refFile.empty()) {
    cerr << "loading reference data:" << endl;
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
