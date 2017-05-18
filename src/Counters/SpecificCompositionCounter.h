#ifndef COUNTERS_SPECIFICCOMPOSITIONCOUNTER_H_
#define COUNTERS_SPECIFICCOMPOSITIONCOUNTER_H_

#include "../Pattern/Pattern.h"
#include <list>
#include <boost/circular_buffer.hpp>
#include <memory>
#include "../Motifs/CompactMotifBuilder.h"

#include "IMotifCounter.h"

class SpecificCompositionCounter: public IMotifCounter {
private:
    std::unique_ptr<unsigned[]> plusAndMinus;
    cell kmer;
    bool rc, fuzzySpacer, fuzzyOrder, fuzzyOrientation;
    unsigned kmersTotal;

    std::vector<Pattern> patterns, rcPatterns;
    std::vector<std::string> strPatterns;

    std::shared_ptr<IDataStructure> result;

    int window, k, end, center, minSpacer, maxSpacer, posOffset;
    CompactMotifBuilder builder;
    boost::circular_buffer<unsigned> buffer;
    boost::circular_buffer<bool> validity;

    void step();
    inline void scroll(unsigned nucleotide);

    inline unsigned getID(unsigned kmer, unsigned pattern_id,
                          int spacer, bool pattern_left) const;
public:
    SpecificCompositionCounter(
        const DataStructureFactory& factory,
        const std::vector<std::string> & geneLabels,
        const std::vector<std::string> & patterns,
        unsigned k,
        int minSpacer,
        int maxSpacer,
        bool rc=false,
        bool fuzzySpacer=false,
        bool fuzzyOrder=false,
        bool fuzzyOrientation=false);

    virtual void initGene(unsigned gene);
    virtual void skip();
    virtual void count(unsigned element);
    virtual void init(const DataStructureFactory& factory,
                      const std::function<std::string (unsigned)> elementLabelGenerator,
                      const std::vector<std::string>& geneLabels);
    virtual void finalizeGene();
    virtual std::shared_ptr<IDataStructure> getResult() const;
    virtual ~SpecificCompositionCounter() {};
};


#endif /* COUNTERS_CLOSETOPATTERNCOUNTER_H_ */
