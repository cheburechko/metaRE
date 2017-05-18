#include "SpecificMotifCounter.h"
#include <algorithm>

SpecificMotifCounter::SpecificMotifCounter(
    const DataStructureFactory& factory,
    const std::vector<std::string> & geneLabels,
    const std::vector<std::string> & patterns,
    bool rc
) :
    rc(rc),
    patterns(),
    pos(0),
    builder(longestPattern(patterns))
{
    std::function<std::string (unsigned)> elementLabelGenerator =
        [patterns](unsigned id){return patterns[id];};
    init(factory, elementLabelGenerator, geneLabels);
    std::for_each(patterns.begin(), patterns.end(),
                   [this](std::string motif){this->patterns.push_back(IUPACMotif(motif));});
}

size_t SpecificMotifCounter::longestPattern(const std::vector<std::string> & patterns) {
    size_t k = 0;
    for (const auto& pattern: patterns) {
        k = std::max(pattern.size(), k);
    }
    return k;
}

void SpecificMotifCounter::initGene(unsigned gene){
    builder.clear();
    pos=0;
    result->sGeneInput(gene);
}

void SpecificMotifCounter::skip() {
    pos++;
    builder.skip();
}

void SpecificMotifCounter::count(unsigned nucleotide){
    builder.putCompact(nucleotide);
    pos++;
    for (unsigned i = 0; i < patterns.size(); i++) {
        if (builder.matches(patterns[i], rc)) {
            result->sElementInput(i, pos);
        }
    }
}

void SpecificMotifCounter::init(const DataStructureFactory& factory,
                                const std::function<std::string (unsigned)> elementLabelGenerator,
                  const std::vector<std::string>& geneLabels) {
    result = std::shared_ptr<IDataStructure>(factory.create(elementLabelGenerator, geneLabels));
}

void SpecificMotifCounter::finalizeGene() {}

std::shared_ptr<IDataStructure> SpecificMotifCounter::getResult() const{
    return result;
}
SpecificMotifCounter::~SpecificMotifCounter() {}
