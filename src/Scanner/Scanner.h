/*
 * MotifCounter.h
 *
 *  Created on: 22.07.2015
 *      Author: thepunchy
 */

#ifndef MOTIFSCANNER_H_
#define MOTIFSCANNER_H_

#include <vector>
#include <string>
#include <map>
#include "../Counters/IMotifCounter.h"

class Scanner {
private:
	IMotifCounter* counter;

public:
	Scanner();
	void setCounter(IMotifCounter* counter);
	const IMotifCounter * getCounter();
	void countMotifs(const std::vector<std::string>& dna,
					 const std::vector<unsigned>& geneIDs,
					 const std::vector<std::string>& geneNames);
	virtual ~Scanner();
};

#endif /* MOTIFSCANNER_H_ */
