#include <Rcpp.h>
#include "unittest.h"
#include "bamutil.h"
#include <time.h>
#include "cluster.h"
namespace GENCORE {
UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= BamUtil::test();
    //passed &= Cluster::test();
    Rcpp::warning("\n==========================\n");
    Rcpp::warning("%s\n\n", passed?"PASSED":"FAILED");
    //printf("\n==========================\n");
    //printf("%s\n\n", passed?"PASSED":"FAILED");
}
}
