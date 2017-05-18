#include <testthat.h>
#include <iostream>

#include "fakeit.hpp"

#include "../Counters/SimpleMotifCounter.h"
#include "../Utils/Utils.h"
#include "../DataStructures/MotifPositionsSparse.h"

using ::fakeit::Verify;
using ::fakeit::Fake;
using ::fakeit::When;
using ::fakeit::Mock;
using ::fakeit::_;

void test_counter(SimpleMotifCounter &counter, unsigned k, unsigned batch) {
    counter.initGene(0);
    unsigned i = 0;
    while (i < batch) {
        i++;
        counter.count(2);
    }
    while (i < batch*2) {
        i++;
        counter.count(1);
    }
    counter.finalizeGene();
}

context("SimpleMotifCounter") {
    test_that("initialization") {
        unsigned k = 4, batch = 5;

        DataStructureFactory factory;

        std::vector<std::string> geneNames({"gene1"});

        test_that("orientation dependent counter works") {
            SimpleMotifCounter counter(factory, geneNames, k, false);
            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            Fake(Method(spy, sElementInput));

            test_counter(counter, k, batch);

            CATCH_CHECK_NOTHROW(Verify(Method(spy, sGeneInput).Using(0)).Exactly(1));
            CATCH_CHECK_NOTHROW(Verify(Method(spy, sElementInput).Using(0x55, _)).Exactly(batch-k+1));
            CATCH_CHECK_NOTHROW(Verify(Method(spy, sElementInput).Using(Utils::reverseComplement(0x55, k), _)).Exactly(batch-k+1));
        };

        test_that("orientation independent counter works")
        {
            SimpleMotifCounter counter(factory, geneNames, k, true);
            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            Fake(Method(spy, sElementInput));

            test_counter(counter, k, batch);

            CATCH_CHECK_NOTHROW(Verify(Method(spy, sGeneInput).Using(0)).Exactly(1));
            CATCH_CHECK_NOTHROW(Verify(Method(spy, sElementInput).Using(0x55, _)).Exactly((batch-k+1)*2));
        };
    };
}
