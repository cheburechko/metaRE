#include "GeneComposition.h"
#include <algorithm>
#include <climits>

GeneComposition::GeneComposition(const std::function<std::string (unsigned)> elementLabelGenerator,
                                           const std::vector<std::string>& geneLabels) :
    elementLabelGenerator(elementLabelGenerator),
    geneLabels(geneLabels),
    curGene(0),
    empty(-1),
    k(INT_MAX)
{}

void GeneComposition::sGeneInput(unsigned gene) {
    curGene = gene;
}

void GeneComposition::sElementInput(unsigned element, int position) {
    if (data[curGene].size() <= position) {
        data[curGene].resize(position+1, empty);
        k = std::min(position, k);
    }
    elements.insert(element);
    data[curGene][position] = element;
}

unsigned GeneComposition::getElementCount() const{
    return elements.size();
}

unsigned GeneComposition::getGeneCount() const{
    return geneLabels.size();
}

std::string GeneComposition::getElementLabel(unsigned element) const{
    return elementLabelGenerator(element);
}

const std::string& GeneComposition::getGeneLabel(unsigned gene) const{
    return geneLabels[gene];
}

std::unordered_map<unsigned, std::string> GeneComposition::getElementLabels() const{
    std::unordered_map<unsigned, std::string> result;
    for (auto it : elements) {
        result[it] = elementLabelGenerator(it);
    }
    return result;
}

const std::vector<std::string>& GeneComposition::getGeneLabels() const{
    return geneLabels;
}

const std::unordered_map<unsigned, std::vector<int>>& GeneComposition::getStructure() const {
    return data;
}

SEXP GeneComposition::getSEXP() const {
    // Available only for oligomer enumeration
    auto motifs = getElementLabels();

    Rcpp::List result;
    for (auto it : data) {
        std::string name = geneLabels[it.first];
        std::vector<Rcpp::String> buf;
        for (auto it2 = it.second.cbegin()+k; it2 != it.second.cend(); it2++) {
            if (*it2 == empty) {
                buf.push_back(NA_STRING);
            }
            buf.push_back(motifs[*it2]);
        }
        Rcpp::CharacterVector sequence(buf.begin(), buf.end());
        result[name] = sequence;
    }
    return result;
}
