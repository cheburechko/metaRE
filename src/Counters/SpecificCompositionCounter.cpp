#include "SpecificCompositionCounter.h"

#include "../Utils/Utils.h"
#include <string>
#include <algorithm>

SpecificCompositionCounter::SpecificCompositionCounter(
    const DataStructureFactory& factory,
    const std::vector<std::string> & geneLabels,
    const std::vector<std::string> & patterns,
    unsigned k,
    int minSpacer,
    int maxSpacer,
    bool rc,
    bool fuzzySpacer,
    bool fuzzyOrder,
    bool fuzzyOrientation
):
    rc(rc),
    result(nullptr),
    minSpacer(minSpacer),
    maxSpacer(maxSpacer),
    fuzzySpacer(fuzzySpacer),
    fuzzyOrder(fuzzyOrder),
    fuzzyOrientation(fuzzyOrientation),
    window(maxSpacer-minSpacer+1),
    k(k),
    builder(k),
    end(0),
    buffer(2*(maxSpacer+k)+1),
    validity(2*(maxSpacer+k)+1),
    center(-maxSpacer-k-1)
{
    unsigned limit = sizeof(unsigned) * 8 / COMPACT_SIZE;
    if (k > limit) {
        throw std::runtime_error("Error: coupling oligomer length has a limit of " + std::to_string(limit));
    }
    plusAndMinus = Utils::getPlusMinusMapping(k);
    kmersTotal = Utils::getKmersTotal(k);

    // Generating labels
    for(unsigned patternID = 0; patternID < patterns.size(); patternID++) {
        std::string strPattern = patterns[patternID];
        this -> patterns.push_back(Pattern(strPattern));
        rcPatterns.push_back(Pattern(Utils::reverseComplement(strPattern)));
    }

    std::function<std::string (unsigned)> elementLabelGenerator =
        [this, patterns](unsigned id) {
            int kmer, pattern, orientation, spacer;
            std::string strSpacer;
            kmer = id % kmersTotal;
            //kmer = this->fuzzyOrientation && this->rc && ? this->plusAndMinus[kmer] : kmer;
            id /= kmersTotal;
            if (this->fuzzySpacer) {
                strSpacer = std::to_string(this->minSpacer) +
                    ".." + std::to_string(this->maxSpacer);
            } else {
                spacer = id % this->window + this->minSpacer;
                id /= this->window;
                strSpacer = std::to_string(spacer);
            }
            pattern = id % patterns.size();
            orientation = id / patterns.size();
            if (!orientation) {
                return patterns[pattern] + "_" + strSpacer +
                    "_" + Utils::intToString(kmer, this->k, true);
            } else {
                return Utils::intToString(kmer, this->k, true) +
                    "_" + strSpacer + "_" + patterns[pattern];
            }
        };
    init(factory, elementLabelGenerator, geneLabels);
}

inline unsigned SpecificCompositionCounter::getID(unsigned kmer,
                                                  unsigned pattern_id,
                                                  int spacer,
                                                  bool pattern_left) const {
    unsigned patternAndOrientation = pattern_id + patterns.size() *
        (fuzzyOrder || pattern_left ? 0 : 1);
    kmer = fuzzyOrientation && rc ? plusAndMinus[kmer] : kmer;
    if (fuzzySpacer) {
        return kmer + kmersTotal * patternAndOrientation;
    }
    return kmer + kmersTotal * (spacer-minSpacer + window * patternAndOrientation);
}

std::shared_ptr<IDataStructure> SpecificCompositionCounter::getResult() const {
    return result;
}

void SpecificCompositionCounter::step() {
    if (!buffer.full()) {
        center++;
        end++;
    }
    builder.write(&kmer);
    buffer.push_back(kmer);
    validity.push_back(builder.ready());
    if (center < 0) {
        return;
    }
    if (validity[center]) {
        for (unsigned patternID = 0; patternID < patterns.size(); patternID++) {
            if (patterns[patternID].check(buffer[center])) {
                for (int i = 0; i < end; i++) {
                    if (i < center && (patterns[patternID].check(buffer[i]) ||
                        (rc && rcPatterns[patternID].check(buffer[i]))
                    )) {
                        continue;
                    }
                    if (abs(center - i) >= k+minSpacer && validity[i]) {
                        result -> sElementInput(
                                getID(
                                    buffer[i],
                                    patternID,
                                    abs(center - i) - k,
                                    i > center
                                ),
                                builder.getPos() - posOffset - end + (i > center ? i : center) + 1
                        );
                    }
                }
            }
            if (rc && (rcPatterns[patternID].check(buffer[center]))) {
                for (int i = 0; i < end; i++) {
                    if (i < center && (patterns[patternID].check(buffer[i]) ||
                        rcPatterns[patternID].check(buffer[i])
                    )) {
                        continue;
                    }
                    if (abs(center - i) >= k+minSpacer && validity[i]) {
                        result -> sElementInput(
                                getID(
                                    Utils::reverseComplement(buffer[i], k),
                                    patternID,
                                    abs(center - i) - k,
                                    i < center
                                ),
                                builder.getPos() - posOffset - end + (i > center ? i : center) + 1);
                    }
                }
            }
        }
    }
}

void SpecificCompositionCounter::initGene(unsigned gene) {
    result -> sGeneInput(gene);
    buffer.clear();
    validity.clear();
    end = 0;
    center = -maxSpacer-k-1;
    builder.clear();
    posOffset = 0;
}

void SpecificCompositionCounter::skip() {
    builder.skip();
    step();
}

void SpecificCompositionCounter::count(unsigned nucleotide) {
    builder.put(nucleotide);
    step();
}

void SpecificCompositionCounter::init(const DataStructureFactory& factory,
                                      const std::function<std::string (unsigned)> elementLabelGenerator,
                                      const std::vector<std::string>& geneLabels) {
    result = std::shared_ptr<IDataStructure>(factory.create(elementLabelGenerator, geneLabels));
}

void SpecificCompositionCounter::finalizeGene() {
    while(end > center) {
        end--;
        posOffset++;
        skip();
    }
}
