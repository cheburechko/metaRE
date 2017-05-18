/*
 * SimpleMotifCounter.h
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#ifndef COUNTERS_SIMPLEMOTIFCOUNTER_H_
#define COUNTERS_SIMPLEMOTIFCOUNTER_H_

#include "../DataStructures/MotifPositionsSparse.h"
#include "IMotifCounter.h"
#include "../Motifs/CompactMotifBuilder.h"
#include <memory>

class SimpleMotifCounter: public IMotifCounter {
private:
	unsigned kmersTotal, k;
    cell kmer;
    bool rc;
	std::unique_ptr<unsigned[]> plusAndMinus;
	std::shared_ptr<IDataStructure> result;
	CompactMotifBuilder builder;
public:
	SimpleMotifCounter(const DataStructureFactory& factory,
                    const std::vector<std::string> & geneLabels,
                       unsigned k, bool rc);
    virtual unsigned getK() const { return k; };
	virtual void initGene(unsigned gene);
	virtual void count(unsigned element);
	virtual void skip();
	virtual void finalizeGene() {};
	virtual void init(const DataStructureFactory& factory,
                   const std::function<std::string (unsigned)> elementLabelGenerator,
                   const std::vector<std::string>& geneLabels);
	virtual std::shared_ptr<IDataStructure> getResult() const;
	virtual ~SimpleMotifCounter() {};
};

#endif /* COUNTERS_SIMPLEMOTIFCOUNTER_H_ */
