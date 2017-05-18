#include "ElementCounts.h"

ElementCounts::ElementCounts(const std::function<std::string (unsigned)> elementLabelGenerator,
                             const std::vector<std::string>& geneLabels) :
    elementLabelGenerator(elementLabelGenerator),
    geneLabels(geneLabels),
    curGene(0)
{}

void ElementCounts::sGeneInput(unsigned gene) {
    curGene = gene;
}

void ElementCounts::sElementInput(unsigned element, int position) {
    if (data.find(element) == data.end()) {
        data[element].push_back(1);
        return;
    }
    data[element][0]++;
}

unsigned ElementCounts::getElementCount() const{
    return data.size();
}

unsigned ElementCounts::getGeneCount() const{
    return geneLabels.size();
}

std::string ElementCounts::getElementLabel(unsigned element) const{
    return elementLabelGenerator(element);
}

const std::string& ElementCounts::getGeneLabel(unsigned gene) const{
    return geneLabels[gene];
}

const std::vector<std::string>& ElementCounts::getGeneLabels() const{
    return geneLabels;
}

const std::unordered_map<unsigned, std::vector<int>>& ElementCounts::getStructure() const {
    return data;
}

SEXP ElementCounts::getSEXP() const {
    std::vector<unsigned> buffer;
    std::vector<std::string> nameBuffer;
    for (auto it : data) {
        buffer.push_back(it.second[0]);
        nameBuffer.push_back(elementLabelGenerator(it.first));
    }
    Rcpp::NumericVector filteredData = Rcpp::wrap(buffer);
    filteredData.attr("names") = nameBuffer;
    return filteredData;
}
