#ifndef MOTIFPOSITIONS_H_
#define MOTIFPOSITIONS_H_

#include "IDataStructure.h"
#include <vector>
#include <map>
#include <string>

class MotifPositions : public IDataStructure {
private:
    const std::vector<std::string> geneLabels;
    const std::function<std::string (unsigned)> elementLabelGenerator;
    unsigned curGene;

    std::unordered_map<unsigned, std::unordered_map<int, std::vector<int> > > data;
    std::unordered_map<unsigned, std::vector<int> > dummy;

public:
    MotifPositions(const std::function<std::string (unsigned)> elementLabelGenerator,
                   const std::vector<std::string> & geneLabels);

    virtual void sGeneInput(unsigned gene);
    virtual void sElementInput(unsigned element, int position);

    virtual const std::unordered_map<unsigned, std::vector<int> >& getStructure() const;
    virtual const std::unordered_map<unsigned, std::unordered_map<int, std::vector<int> > >& getPositions() const;
    virtual SEXP getSEXP() const;

    virtual unsigned getElementCount() const;
    virtual unsigned getGeneCount() const;
    virtual std::string getElementLabel(unsigned element) const;
    virtual const std::string& getGeneLabel(unsigned gene) const;
    virtual std::unordered_map<unsigned, std::string> getElementLabels() const;
    virtual const std::vector<std::string> & getGeneLabels() const;
    virtual ~MotifPositions() {};
};

#endif /* MOTIFPOSITIONS_H_ */
