#include "Scanner.h"

#include <list>
#include <cmath>
#include <algorithm>
#include "../Utils/Utils.h"

using namespace std;

Scanner::Scanner() {}

void Scanner::countMotifs(const std::vector<std::string>& dna,
		const std::vector<unsigned>& geneIDs,
		const std::vector<std::string>& geneNames) {
	for (unsigned geneNumber = 0; geneNumber < dna.size(); geneNumber++) {
		counter->initGene(geneIDs[geneNumber]);

		for (unsigned i = 0; i < dna[geneNumber].size(); i++) {
		    unsigned nucleotide = Utils::charToInt(dna[geneNumber][i]);
			if (nucleotide != -1) {
				counter->count(nucleotide);
			} else {
				counter->skip();
			}
		}
		counter->finalizeGene();
	}
}

void Scanner::setCounter(IMotifCounter* counter) {
	this->counter = counter;
}

const IMotifCounter* Scanner::getCounter() {
	return counter;
}

Scanner::~Scanner() {
}
