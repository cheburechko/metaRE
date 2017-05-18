#ifndef COUNTERS_REPEATCOUNTER_H_
#define COUNTERS_REPEATCOUNTER_H_

#include "../DataStructures/DataStructureFactory.h"
#include "IMotifCounter.h"
#include "../Motifs/CompactMotifBuilder.h"
#include "MotifBuffer.hpp"
#include <memory>

class RepeatCounter: public IMotifCounter {
private:
    int kmersTotal, k, window, minSpacer, maxSpacer;
    CompactMotifBuilder builder;
    MotifBuffer<cell> buffer;
    std::vector<unsigned> presentKmerMap;

    std::unique_ptr<unsigned[]> plusAndMinus;
    std::shared_ptr<IDataStructure> result;

    void step();

public:
    unsigned getID(unsigned kmer, unsigned spacer, unsigned orientation);

    static const unsigned DIRECT = 0, INVERTED = 1, EVERTED = 2;
    RepeatCounter(
        const DataStructureFactory& factory,
        const std::vector<std::string> & geneLabels,
        unsigned k,
        int minSpacer,
        int maxSpacer
    );
    virtual void initGene(unsigned gene);
    virtual void count(unsigned nucleotide);
    virtual void skip();
    virtual void finalizeGene();
    virtual void init(const DataStructureFactory& factory,
                      const std::function<std::string (unsigned)> elementLabelGenerator,
                      const std::vector<std::string>& geneLabels);
    virtual std::shared_ptr<IDataStructure> getResult() const;
    virtual ~RepeatCounter() {};
};

#endif /* COUNTERS_REPEATCOUNTER_H_ */
