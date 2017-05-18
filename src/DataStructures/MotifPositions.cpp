#include "MotifPositions.h"

MotifPositions::MotifPositions(const std::function<std::string (unsigned)> elementLabelGenerator,
                               const std::vector<std::string>& geneLabels) :
    elementLabelGenerator(elementLabelGenerator),
    geneLabels(geneLabels),
    curGene(0)
{}

void MotifPositions::sGeneInput(unsigned gene) {
    curGene = gene;
}

void MotifPositions::sElementInput(unsigned element, int position) {
    auto it = data[element].find(curGene);
    if (it != data[element].end()) {
        (it->second).push_back(position);
    } else {
        data[element].insert(
            std::pair<int, std::vector<int> >(curGene, std::vector<int>({position}))
        );
    }
}

unsigned MotifPositions::getElementCount() const{
    return data.size();
}

unsigned MotifPositions::getGeneCount() const{
    return geneLabels.size();
}

std::string MotifPositions::getElementLabel(unsigned element) const{
    return elementLabelGenerator(element);
}

std::unordered_map<unsigned, std::string> MotifPositions::getElementLabels() const {
    std::unordered_map<unsigned, std::string> result;
    for (auto it : data) {
        result[it.first] = getElementLabel(it.first);
    }
    return result;
}


const std::string& MotifPositions::getGeneLabel(unsigned gene) const{
    return geneLabels[gene];
}

const std::vector<std::string>& MotifPositions::getGeneLabels() const{
    return geneLabels;
}

const std::unordered_map<unsigned, std::vector<int> >& MotifPositions::getStructure() const {
    return dummy;
}

const std::unordered_map<unsigned, std::unordered_map<int, std::vector<int> > >& MotifPositions::getPositions() const {
    return data;
}

SEXP MotifPositions::getSEXP() const {
    Rcpp::List result;

    for (auto it : data) {
        Rcpp::List sublist;
        for (auto it2 : it.second) {
            sublist[geneLabels[it2.first]] = Rcpp::wrap(it2.second);
        }
        std::string element = elementLabelGenerator(it.first);
        result[element] = sublist;
    }
    result.attr("genes") = geneLabels;
    return result;
}
