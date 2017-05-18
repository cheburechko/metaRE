/*
 * SimpleMotifCounter.cpp
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#include "SimpleMotifCounter.h"
#include "../Utils/Utils.h"

SimpleMotifCounter::SimpleMotifCounter(
    const DataStructureFactory& factory,
    const std::vector<std::string> & geneLabels,
    unsigned k,
    bool rc
) :
    rc(rc),
	result(nullptr),
	k(k),
	builder(k)
{
    unsigned limit = sizeof(unsigned) * 8 / COMPACT_SIZE;
    if (k > limit) {
        throw std::runtime_error("Error: oligomer length has a limit of " + std::to_string(limit));
    }
	plusAndMinus = Utils::getPlusMinusMapping(k);
	kmersTotal = Utils::getKmersTotal(k);
	std::function<std::string (unsigned)> elementLabelGenerator =
	    [this](unsigned id) {
	        std::string label = Utils::intToString(id, this->k, true);
	        if (this->rc) {
	            label += " | " + Utils::reverseComplement(label, true);
	        }
	        return label;
	    };
	init(factory, elementLabelGenerator, geneLabels);
}

void SimpleMotifCounter::count(unsigned nucleotide) {
    builder.put(nucleotide);
	if (builder.ready()) {
	    builder.write(&kmer);
	    if (rc) {
	        result->sElementInput(plusAndMinus[kmer], builder.getPos());
	    } else {
	        result->sElementInput(kmer, builder.getPos());
	    }
	}
}

void SimpleMotifCounter::skip() {
    builder.skip();
}

std::shared_ptr<IDataStructure> SimpleMotifCounter::getResult() const {
	return result;
}

void SimpleMotifCounter::init(const DataStructureFactory& factory,
                              const std::function<std::string (unsigned)> elementLabelGenerator,
                              const std::vector<std::string>& geneLabels) {
    result = std::shared_ptr<IDataStructure>(factory.create(elementLabelGenerator, geneLabels));
}

void SimpleMotifCounter::initGene(unsigned gene) {
    builder.clear();
	result->sGeneInput(gene);
}
