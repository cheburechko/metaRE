#include <Rcpp.h>
#include "stat_tests.hpp"
using namespace Rcpp;

inline unsigned getIndex(int element, int experiment, int elements) {
    return element + experiment*elements;
}

// [[Rcpp::export]]
NumericMatrix massFisherTest(const LogicalMatrix& experiments, const IntegerVector& sums,
                             const List& elements, std::string altString) {
    // Initializing
    Alternative alternative = strToAlternative(altString);
    unsigned resultSize = elements.size() * experiments.ncol();
    NumericMatrix result(elements.size(), experiments.ncol());

    unsigned totalGenes = experiments.nrow();
    unsigned elem, noElem, elemUp, noElemUp;

    for (unsigned element = 0; element < elements.size(); element++ ) {
        int regGenes = LENGTH(elements[element]);
        std::vector<int> geneIndices = as<std::vector<int> >(elements[element]);
        for (unsigned experiment = 0; experiment < experiments.ncol(); experiment++) {
            unsigned index = getIndex(element, experiment, elements.size());
            elem = regGenes;
            noElem = totalGenes - regGenes;
            elemUp = 0;
            for (unsigned gene = 0; gene < regGenes; gene++) {
                if (experiments(geneIndices[gene]-1, experiment)) {
                    elemUp++;
                }
            }
            noElemUp = sums[experiment] - elemUp;
            result[index] = fisherTest(elemUp, elem, noElemUp, noElem, alternative);
        }
    }
    return result;
}


// [[Rcpp::export]]
NumericVector quickFisherTest(NumericVector eff1, NumericVector n1, NumericVector eff2, NumericVector n2,
                            std::string alternative){
    NumericVector result(eff1.length());
    for (unsigned i = 0; i < eff1.length(); i++) {
        result[i] = fisherTest(eff1[i], n1[i], eff2[i], n2[i], strToAlternative(alternative));
    }
    return result;
}
