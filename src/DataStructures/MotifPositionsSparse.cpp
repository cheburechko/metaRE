/*
 * MotifPositions.cpp
 *
 *  Created on: 22.07.2015
 *      Author: thepunchy
 */

#include "MotifPositionsSparse.h"

MotifPositionsSparse::MotifPositionsSparse(const std::function<std::string (unsigned)> elementLabelGenerator,
										   const std::vector<std::string>& geneLabels, Rcpp::Function createGCS) :
		elementLabelGenerator(elementLabelGenerator),
		geneLabels(geneLabels),
		curGene(0),
		createGCS(createGCS)
{}

void MotifPositionsSparse::sGeneInput(unsigned gene) {
	curGene = gene;
}

void MotifPositionsSparse::sElementInput(unsigned element, int position) {
    if (data[element].empty() || data[element].back() != curGene) {
	    data[element].push_back(curGene);
    }
}

unsigned MotifPositionsSparse::getElementCount() const{
	return data.size();
}

unsigned MotifPositionsSparse::getGeneCount() const{
	return geneLabels.size();
}

std::string MotifPositionsSparse::getElementLabel(unsigned element) const{
	return elementLabelGenerator(element);
}

const std::string& MotifPositionsSparse::getGeneLabel(unsigned gene) const{
	return geneLabels[gene];
}

const std::vector<std::string>& MotifPositionsSparse::getGeneLabels() const{
	return geneLabels;
}

const std::unordered_map<unsigned, std::vector<int>>& MotifPositionsSparse::getStructure() const {
    return data;
}

SEXP MotifPositionsSparse::getSEXP() const {
    Rcpp::CharacterVector genes = Rcpp::wrap(geneLabels);

    std::unordered_map<std::string, std::vector<int>> cppResult;
    for (auto it : data) {
        std::string elementName = elementLabelGenerator(it.first);
        cppResult[elementName] = std::vector<int>(it.second.size());
        std::transform(
            it.second.begin(), it.second.end(),
            cppResult[elementName].begin(),
            [](int x){return x+1;}
        );
    }

    Rcpp::List result = Rcpp::wrap(cppResult);
    return createGCS(result, genes);
}
