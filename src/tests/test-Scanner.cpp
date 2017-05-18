#include <testthat.h>
#include <iostream>
#include "fakeit.hpp"

#include "../Counters/IMotifCounter.h"
#include "../Scanner/Scanner.h"
#include "../Utils/Utils.h"

using ::fakeit::Verify;
using ::fakeit::VerifyNoOtherInvocations;
using ::fakeit::Fake;
using ::fakeit::Mock;
using ::fakeit::_;

context("Scanner") {
    test_that("initialize") {
        Mock<IMotifCounter> counter;
        Scanner scanner;
        std::vector<std::string> dnas({
            "aaaaaaaaaa", //10
            "cccccccccccccc", //14
            "ttttttttxttttttttt" //18
        });
        std::vector<std::string> geneNames({"gene1", "gene2", "gene3"});
        std::vector<unsigned> geneIDs({0,1,2});


        Fake(Method((counter), count));
        Fake(Method((counter), initGene));
        Fake(Method((counter), finalizeGene));
        Fake(Method((counter), skip));
        scanner.setCounter(&(counter.get()));

        test_that("can set and get counter")
        {
            auto result = scanner.getCounter();
            expect_true(result == &(counter.get()));
        }

        test_that("scan with no parameters works")
        {
            scanner.countMotifs(dnas, geneIDs, geneNames);

            Verify(Method((counter), count).Using(Utils::charToInt('a')))
                .Exactly(dnas[0].length());
            Verify(Method((counter), count).Using(Utils::charToInt('c')))
                .Exactly(dnas[1].length());
            Verify(Method((counter), count).Using(Utils::charToInt('t')))
                .Exactly(dnas[2].length()-1);
            Verify(Method((counter), skip))
                .Exactly(1);
        }
    }
}
