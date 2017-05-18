#include "DataStructureFactory.h"
#include "ElementCounts.h"
#include "MotifPositions.h"
#include "MotifPositionsSparse.h"
#include "GeneComposition.h"
#include <stdexcept>

const std::map<std::string, DataStructureFactory::type> DataStructureFactory::dataTypeDict = {
    {"counts", DataStructureFactory::type::ElementCounts},
    {"genes", DataStructureFactory::type::MotifPositionsSparse},
    {"positions", DataStructureFactory::type::MotifPositions},
    {"composition", DataStructureFactory::type::GeneComposition}
};

DataStructureFactory::DataStructureFactory() :
    _type(DataStructureFactory::type::MotifPositionsSparse),
    createGCS("list")
{}

DataStructureFactory::DataStructureFactory(const DataStructureFactory& other) :
    _type(other._type),
    createGCS(other.createGCS)
{}

void DataStructureFactory::setType(DataStructureFactory::type value) {
    _type = value;
}

bool DataStructureFactory::setType(std::string value) {
    auto it = dataTypeDict.find(value);
    if (it == dataTypeDict.end()) {
        return false;
    }
    _type = it -> second;
    return true;
}

DataStructureFactory::type DataStructureFactory::getType() const {
    return _type;
}

void DataStructureFactory::setCreateGCS(Rcpp::Function func) {
    createGCS = func;
}

IDataStructure * DataStructureFactory::create(
        const std::function<std::string (unsigned)> elementLabelGenerator,
        const std::vector<std::string>& geneLabels) const {
    switch(_type) {
        case DataStructureFactory::type::MotifPositionsSparse:
            return new MotifPositionsSparse(elementLabelGenerator, geneLabels, createGCS);
        case DataStructureFactory::type::ElementCounts:
            return new ElementCounts(elementLabelGenerator, geneLabels);
        case DataStructureFactory::type::MotifPositions:
            return new MotifPositions(elementLabelGenerator, geneLabels);
        case DataStructureFactory::type::GeneComposition:
            return new GeneComposition(elementLabelGenerator, geneLabels);
        default:
            throw std::invalid_argument("Illegal data structure type");
    }
}
