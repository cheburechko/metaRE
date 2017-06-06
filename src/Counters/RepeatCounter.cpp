#include "RepeatCounter.h"
#include "../Utils/Utils.h"

RepeatCounter::RepeatCounter(
    const DataStructureFactory& factory,
    const std::vector<std::string> & geneLabels,
    unsigned k,
    int minSpacer,
    int maxSpacer
) :
    result(nullptr),
    k(k),
    minSpacer(minSpacer),
    maxSpacer(maxSpacer),
    window(maxSpacer-minSpacer),
    builder(k),
    buffer(2*(maxSpacer+k)+1)
{
    plusAndMinus = Utils::getPlusMinusMapping(k);
    kmersTotal = Utils::getKmersTotal(k);
    std::vector<std::string> elementLabels = Utils::generateAllKmerLabels(k, true);

    presentKmerMap = std::vector<unsigned>(kmersTotal, 0);

    const std::function<std::string (unsigned)> elementLabelGenerator =
        [this](unsigned id){
            unsigned kmer = id % this->kmersTotal;
            int spacer = this->minSpacer + (int)(id / this->kmersTotal % (this->window+1));
            unsigned orientation = id / this->kmersTotal / (this->window+1);
            std::string label = Utils::intToString(kmer, this->k, true);
            switch(orientation) {
            case DIRECT:
                return label + '_' + std::to_string(spacer) + '_' + label;
            case INVERTED:
                return Utils::reverseComplement(label,true) + '_' + std::to_string(spacer) + '_' + label;
            case EVERTED:
                return label + '_' + std::to_string(spacer) + '_' + Utils::reverseComplement(label, true);
            }
            return std::string("NULL");
        };

    init(factory, elementLabelGenerator, geneLabels);
}

unsigned RepeatCounter::getID(unsigned kmer, unsigned spacer,
                              unsigned orientation)
{
    return plusAndMinus[kmer] + kmersTotal*(spacer-minSpacer+orientation*(window+1));
};

void RepeatCounter::step() {
    cell kmer;
    bool valid;
    if (builder.ready()) {
        kmer = 0;
        builder.write(&kmer);
        presentKmerMap[plusAndMinus[kmer]]++;

        valid = buffer.isValid(0);
        buffer.put(kmer);
    } else {
        valid = buffer.isValid(0);
        buffer.skip();
    }

    cell center = buffer.getCenter();
    if (buffer.centerAvailable()) {
        presentKmerMap[plusAndMinus[center]]--;
    }

    if (buffer.centerAvailable() && presentKmerMap[plusAndMinus[center]] > 0) {
        for (int i = buffer.getCenterPos()+k+minSpacer; i < buffer.size(); i++) {
            int spacer = buffer.posFromCenter(i) - k;
            if (buffer.isValid(i)) {
                cell kmer = buffer.get(i);
                if (kmer == center) {
                    result->sElementInput(
                            getID(kmer, spacer, DIRECT),
                            buffer.absPosition(i+1)
                    );
                } else if (kmer == plusAndMinus[center] || plusAndMinus[kmer] == center) {
                    bool orig = center == plusAndMinus[center];
                    result->sElementInput(
                            getID(kmer, spacer, orig ? EVERTED : INVERTED),
                            buffer.absPosition(i+1)
                    );
                }
            }
        }
    }
}

void RepeatCounter::count(unsigned nucleotide) {
    builder.put(nucleotide);
    step();


}

void RepeatCounter::skip() {
    builder.skip();
    step();
}

std::shared_ptr<IDataStructure> RepeatCounter::getResult() const {
    return result;
}

void RepeatCounter::init(const DataStructureFactory& factory,
                         const std::function<std::string (unsigned)> elementLabelGenerator,
                         const std::vector<std::string>& geneLabels) {
    result = std::shared_ptr<IDataStructure>(factory.create(elementLabelGenerator, geneLabels));
}

void RepeatCounter::initGene(unsigned gene) {
    buffer.reset();
    builder.clear();
    result->sGeneInput(gene);
}

void RepeatCounter::finalizeGene() {
    for (int i = 0; i < buffer.size(); i++) {
        skip();
    }
};
