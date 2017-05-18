#include <testthat.h>
#include <iostream>

#include "fakeit.hpp"

#include "../Counters/RepeatCounter.h"
#include "../Utils/Utils.h"
#include "../DataStructures/MotifPositionsSparse.h"

using ::fakeit::Verify;
using ::fakeit::VerifyNoOtherInvocations;
using ::fakeit::Fake;
using ::fakeit::When;
using ::fakeit::Mock;
using ::fakeit::_;

void put_kmer(RepeatCounter &counter,
              unsigned k, unsigned kmer) {
    for (unsigned i = 0; i < k; i++) {
        counter.count((kmer >> ((k-i-1)*2)) & 3);
    }
}

void skip(RepeatCounter &counter, unsigned window) {
    for (unsigned i = 0; i < window; i++) {
        counter.count(0);
    }
}

context("RepeatCounter") {
    test_that("initialize") {
        unsigned k = 4;
        unsigned window = 4;
        unsigned kmersTotal = Utils::getKmersTotal(k);
        std::vector<std::string> geneNames({"gene1", "gene2", "gene3"});

        DataStructureFactory factory;

        test_that("counter counts fine") {
            RepeatCounter counter(factory, geneNames, k, 0, window);
            Mock<IDataStructure> data(*(counter.getResult()));
            Fake(Method(data, sGeneInput));
            When(Method(data, sElementInput)).AlwaysDo(
                [](unsigned id, int pos){}
            );

            counter.initGene(0);
            put_kmer(counter, k, Utils::stringToInt("aaaa"));
            put_kmer(counter, k, Utils::stringToInt("cccc"));
            put_kmer(counter, k, Utils::stringToInt("gggg"));
            put_kmer(counter, k, Utils::stringToInt("tttt"));
            put_kmer(counter, k, Utils::stringToInt("aaaa"));
            put_kmer(counter, k, Utils::stringToInt("cccc"));
            put_kmer(counter, k, Utils::stringToInt("gggg"));
            put_kmer(counter, k, Utils::stringToInt("tttt"));
            counter.finalizeGene();

            counter.initGene(1);
            put_kmer(counter, k, Utils::stringToInt("aaaa"));
            put_kmer(counter, k, Utils::stringToInt("tttt"));
            put_kmer(counter, k, Utils::stringToInt("aaaa"));
            put_kmer(counter, k, Utils::stringToInt("tttt"));
            counter.finalizeGene();

            CATCH_CHECK_NOTHROW(Verify(
                Method(data, sGeneInput).Using(0) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aacc"), 4, RepeatCounter::INVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("accc"), 2, RepeatCounter::INVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("cccc"), 0, RepeatCounter::INVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aacc"), 4, RepeatCounter::EVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaac"), 2, RepeatCounter::EVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaaa"), 0, RepeatCounter::EVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aacc"), 4, RepeatCounter::INVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("accc"), 2, RepeatCounter::INVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("cccc"), 0, RepeatCounter::INVERTED), _) +

                Method(data, sGeneInput).Using(1) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaaa"), 0, RepeatCounter::INVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaaa"), 4, RepeatCounter::DIRECT), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaat"), 4, RepeatCounter::DIRECT), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aatt"), 4, RepeatCounter::DIRECT), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaat"), 2, RepeatCounter::EVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("attt"), 4, RepeatCounter::DIRECT), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaaa"), 0, RepeatCounter::EVERTED), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaaa"), 4, RepeatCounter::DIRECT), _) +
                Method(data, sElementInput).Using(counter.getID(Utils::stringToInt("aaaa"), 0, RepeatCounter::INVERTED), _)

            ));

            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(data));
        }
    }
}
