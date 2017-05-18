/*
 * AbstractCounter.h
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#ifndef COUNTERS_IMOTIFCOUNTER_H_
#define COUNTERS_IMOTIFCOUNTER_H_

#include "../DataStructures/DataStructureFactory.h"
#include <memory>
#include <map>

class IMotifCounter {
public:

	IMotifCounter() {};
	virtual void initGene(unsigned gene) = 0;
	virtual void count(unsigned nucleotide) = 0;
	virtual void skip() = 0;
	virtual void finalizeGene() = 0;
	virtual void init(const DataStructureFactory& factory,
	                  const std::function<std::string (unsigned)> elementLabelGenerator,
                      const std::vector<std::string>& geneLabels) = 0;
	virtual std::shared_ptr<IDataStructure> getResult() const = 0;
	virtual ~IMotifCounter() {};
};

#endif /* COUNTERS_IMOTIFCOUNTER_H_ */
