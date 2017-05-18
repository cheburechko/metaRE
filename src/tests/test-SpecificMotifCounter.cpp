#include <testthat.h>
#include <iostream>

#include "fakeit.hpp"

#include "../Counters/SpecificMotifCounter.h"
#include "../Utils/Utils.h"
#include "../DataStructures/MotifPositionsSparse.h"

using ::fakeit::Verify;
using ::fakeit::VerifyNoOtherInvocations;
using ::fakeit::Fake;
using ::fakeit::When;
using ::fakeit::Mock;
using ::fakeit::_;

void test_counter(IMotifCounter &counter) {
    counter.initGene(0);
    counter.count(CHAR_TO_COMPACT('a'));
    counter.count(CHAR_TO_COMPACT('a'));
    counter.count(CHAR_TO_COMPACT('c'));
    counter.count(CHAR_TO_COMPACT('c'));
    counter.count(CHAR_TO_COMPACT('g'));
    counter.count(CHAR_TO_COMPACT('c'));
    counter.count(CHAR_TO_COMPACT('t'));
    counter.count(CHAR_TO_COMPACT('a'));
    counter.count(CHAR_TO_COMPACT('c'));
    counter.finalizeGene();
    counter.initGene(1);
    counter.count(CHAR_TO_COMPACT('g'));
    counter.count(CHAR_TO_COMPACT('t'));
    counter.count(CHAR_TO_COMPACT('a'));
    counter.count(CHAR_TO_COMPACT('g'));
    counter.count(CHAR_TO_COMPACT('c'));
    counter.count(CHAR_TO_COMPACT('g'));
    counter.count(CHAR_TO_COMPACT('g'));
    counter.count(CHAR_TO_COMPACT('t'));
    counter.count(CHAR_TO_COMPACT('t'));
    counter.finalizeGene();
}

context("SpecificMotifCounter") {
    test_that("initialization") {
        DataStructureFactory factory;

        std::vector<std::string> geneNames({"gene1"});

        std::vector<std::string> patterns({"aacc", "rrss", "gtagc"});

        test_that("counter with no flags works") {
            SpecificMotifCounter counter(factory, geneNames, patterns, false);
            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            Fake(Method(spy, sElementInput));

            test_counter(counter);

            CATCH_CHECK_NOTHROW(Verify(
                        Method(spy, sGeneInput).Using(0) +
                        Method(spy, sElementInput).Using(0, _) +
                        Method(spy, sElementInput).Using(1, _) +
                        Method(spy, sGeneInput).Using(1) +
                        Method(spy, sElementInput).Using(2, _) +
                        Method(spy, sElementInput).Using(1, _)
            ));
            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(Method(spy, sElementInput)));
        };

        test_that("counter with pm flag works")
        {
            SpecificMotifCounter counter(factory, geneNames, patterns, true);
            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            Fake(Method(spy, sElementInput));

            test_counter(counter);

            CATCH_CHECK_NOTHROW(Verify(
                        Method(spy, sGeneInput).Using(0) +
                        Method(spy, sElementInput).Using(0, _) +
                        Method(spy, sElementInput).Using(1, _) * 2 +
                        Method(spy, sElementInput).Using(2, _) +
                        Method(spy, sGeneInput).Using(1) +
                        Method(spy, sElementInput).Using(2, _) +
                        Method(spy, sElementInput).Using(1, _) +
                        Method(spy, sElementInput).Using(0, _) +
                        Method(spy, sElementInput).Using(1, _)
            ));
        };
    };
}
