#ifndef COUNTERS_SPECIFICMOTIFCOUNTER_H_
#define COUNTERS_SPECIFICMOTIFCOUNTER_H_

#include "../Motifs/IUPACMotifBuilder.h"
#include "IMotifCounter.h"

class SpecificMotifCounter: public IMotifCounter {
private:
    std::unique_ptr<unsigned[]> plusAndMinus;
    unsigned pos;
    bool rc;
    IUPACMotifBuilder builder;

    std::vector<IUPACMotif> patterns;

    std::shared_ptr<IDataStructure> result;

    static size_t longestPattern(const std::vector<std::string> & patterns);
public:

    SpecificMotifCounter(
        const DataStructureFactory& factory,
        const std::vector<std::string> & geneLabels,
        const std::vector<std::string> & patterns,
        bool rc
    );

    virtual void initGene(unsigned gene);
    virtual void skip();
    virtual void count(unsigned element);
    virtual void init(const DataStructureFactory& factory,
                      const std::function<std::string (unsigned)> elementLabelGenerator,
                      const std::vector<std::string>& geneLabels);
    virtual void finalizeGene();
    virtual std::shared_ptr<IDataStructure> getResult() const;
    virtual ~SpecificMotifCounter();
};


#endif /* COUNTERS_SPECIFICMOTIFCOUNTER_H_ */
