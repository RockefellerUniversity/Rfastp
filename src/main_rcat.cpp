#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#define BUF_SIZE (1<<20)
using namespace Rcpp;

void appendFile(std::string const& outFile, std::string const& inFile) {
    std::ofstream out(outFile, std::ios_base::app |
                               std::ios_base::binary |
                               std::ios_base::out);

    std::ifstream in(inFile, std::ios_base::binary |
                             std::ios_base::in);

    std::vector<char> buffer(BUF_SIZE);
    while (in.read(&buffer[0], buffer.size())) {
        out.write(&buffer[0], buffer.size());
    }

    // Fails when "read" encounters EOF,
    // but potentially still writes *some* bytes to buffer!
    out.write(&buffer[0], in.gcount());
}

// [[Rcpp::export]]
int rcat(std::string output="", Rcpp::List inputFiles = List::create(), int numInFile = 0) {
    if (numInFile < 2) {
        Rcpp::stop("Please specify more than 1 input file.");
    }

    for (int i=0; i<numInFile; i++) {
        appendFile(output, inputFiles[i]);
    }
}
