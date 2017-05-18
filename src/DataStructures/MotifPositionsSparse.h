/*
 * MotifPositions.h
 *
 *  Created on: 22.07.2015
 *      Author: thepunchy
 */

#ifndef MOTIFPOSITIONSSPARSE_H_
#define MOTIFPOSITIONSSPARSE_H_

#include "IDataStructure.h"
#include <Rcpp.h>
#include <vector>
#include <string>
#include <functional>

class MotifPositionsSparse : public IDataStructure {
private:
	const std::vector<std::string> geneLabels;
    const std::function<std::string (unsigned)> elementLabelGenerator;
	unsigned curGene;
	Rcpp::Function createGCS;

	std::unordered_map<unsigned, std::vector<int>> data;

public:
	MotifPositionsSparse(const std::function<std::string (unsigned)> elementLabelGenerator,
						 const std::vector<std::string> & geneLabels,
						 Rcpp::Function createGCS=Rcpp::Function("list"));

    virtual void sGeneInput(unsigned gene);
    virtual void sElementInput(unsigned element, int position);

	virtual const std::unordered_map<unsigned, std::vector<int>>& getStructure() const;
	virtual SEXP getSEXP() const;

	virtual unsigned getElementCount() const;
	virtual unsigned getGeneCount() const;
	virtual std::string getElementLabel(unsigned element) const;
	virtual const std::string& getGeneLabel(unsigned gene) const;
	virtual const std::vector<std::string> & getGeneLabels() const;
	virtual ~MotifPositionsSparse() {};
};

#endif /* MOTIFPOSITIONSSPARSE_H_ */
