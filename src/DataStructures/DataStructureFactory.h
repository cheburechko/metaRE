#ifndef DATASTRUCTUREFACTORY_H_
#define DATASTRUCTUREFACTORY_H_

#include "IDataStructure.h"
#include <functional>
#include <map>

class DataStructureFactory {
public:
    enum class type {ElementCounts, MotifPositionsSparse, MotifPositions, GeneComposition};
    static const std::map<std::string, DataStructureFactory::type> dataTypeDict;
    DataStructureFactory();
    DataStructureFactory(const DataStructureFactory&);

    void setType(type);
    bool setType(std::string);
    type getType() const;

    void setCreateGCS(Rcpp::Function func);

    virtual IDataStructure * create(const std::function<std::string (unsigned)> elementLabelGenerator,
                            const std::vector<std::string>& geneLabels) const;
    virtual ~DataStructureFactory() {};

private:
    type _type;
    Rcpp::Function createGCS;
};

#endif /* DATASTRUCTUREFACTORY_H_ */
