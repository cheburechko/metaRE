#ifndef COUNTERS_MOTIFCOUNTERFACTORY_H_
#define COUNTERS_MOTIFCOUNTERFACTORY_H_

#include <Rcpp.h>
#include "Counters.h"

class MotifCounterFactory {
public:
    MotifCounterFactory() {};
    ~MotifCounterFactory() {};
    IMotifCounter* getMotifCounter(
            const Rcpp::List& counterParams,
            const DataStructureFactory& factory,
            const std::vector<std::string> & geneNames
    ) {
        std::string mode = counterParams["mode"];
        if (!mode.compare("spaced_dyad")) {
            auto pattern = counterParams["patterns"];
            unsigned maxSpacer = counterParams["maxSpacer"];
            unsigned minSpacer = counterParams["minSpacer"];
            bool rc = counterParams["rc"];
            unsigned k = counterParams["k"];
            bool fuzzySpacer = counterParams["fuzzySpacer"];
            bool fuzzyOrder = counterParams["fuzzyOrder"];
            bool fuzzyOrientation = counterParams["fuzzyOrientation"];

            return new SpecificCompositionCounter(
                    factory, geneNames, pattern, k, minSpacer, maxSpacer, rc,
                    fuzzySpacer, fuzzyOrder, fuzzyOrientation
            );
        } else if (!mode.compare("simple")) {
            bool rc = counterParams["rc"];
            unsigned k = counterParams["k"];

            return new SimpleMotifCounter(factory, geneNames, k, rc);
        } else if (!mode.compare("specific_single")) {
            auto pattern = counterParams["patterns"];
            bool rc = counterParams["rc"];

            return new SpecificMotifCounter(factory, geneNames, pattern, rc);
        } else if (!mode.compare("repeat")) {
            unsigned minSpacer = counterParams["minSpacer"];
            unsigned maxSpacer = counterParams["maxSpacer"];
            unsigned k = counterParams["k"];
            return new RepeatCounter(factory, geneNames, k, minSpacer, maxSpacer);
        } else {
            return nullptr;
        }
    }
};

#endif /* COUNTERS_MOTIFCOUNTERFACTORY_H_ */
