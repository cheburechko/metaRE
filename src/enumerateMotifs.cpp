#include <vector>
#include <map>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <Rcpp.h>

#include "Counters/MotifCounterFactory.hpp"
#include "Pattern/Pattern.h"
#include "Scanner/Scanner.h"
#include "DataStructures/MotifPositions.h"
using namespace Rcpp;
using namespace std;

//' @useDynLib metaRE
//' @import Rcpp
// [[Rcpp::export]]
SEXP enumerateMotifsCpp(List parameters, Function createGCS, Function logDebug) {
    map<std::string, unsigned> geneMap;
    vector<unsigned> geneIDs;
    vector<std::string> genes, geneNames;

    logDebug("Initializing parameters...");
    CharacterVector rGenes = parameters["regulatoryRegions"];
    genes = as<std::vector<std::string>>(rGenes);
    if (parameters.containsElementNamed("geneNames")) {
        Vector<STRSXP> rGeneNames = parameters["geneNames"];
        Vector<STRSXP> levels = sort_unique(rGeneNames);

        geneNames = as<std::vector<std::string>>(levels);
        geneIDs = as<std::vector<unsigned>>(match(rGeneNames, levels));
    } else {
        IntegerVector IDs = seq_len(genes.size()) - 1;

        geneNames = as<std::vector<std::string>>(rGenes.names());
        geneIDs = as<std::vector<unsigned>>(IDs);
    }

    DataStructureFactory factory;
    MotifCounterFactory motifCounterFactory;
    string dataType = parameters["data"];
    if (!factory.setType(dataType)) {
        Rprintf("Data type does not exist: %s", dataType.c_str());
        return R_NilValue;
    }
    factory.setCreateGCS(createGCS);

    List counterParams = parameters["counter"];
    auto counter = std::unique_ptr<IMotifCounter>(
        motifCounterFactory.getMotifCounter(counterParams, factory, geneNames)
    );
    if (counter.get() == nullptr) {
        Rprintf("Mode does not exist");
        return R_NilValue;
    }

    logDebug("Done. Computing...");
    Scanner scanner;
    scanner.setCounter(counter.get());
    scanner.countMotifs(genes, geneIDs, geneNames);
    logDebug("Done. Creating structure...");

    SEXP result = counter->getResult()->getSEXP();
    logDebug("Finished.");
    return result;
}
