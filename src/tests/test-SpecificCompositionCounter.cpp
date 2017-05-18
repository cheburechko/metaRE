#include <testthat.h>
#include <iostream>

#include "fakeit.hpp"

#include "../Counters/SpecificCompositionCounter.h"
#include "../Utils/Utils.h"
#include "../DataStructures/MotifPositionsSparse.h"

using ::fakeit::Verify;
using ::fakeit::VerifyNoOtherInvocations;
using ::fakeit::Fake;
using ::fakeit::When;
using ::fakeit::Mock;
using ::fakeit::_;

void put_kmer(SpecificCompositionCounter &counter, unsigned k, unsigned kmer) {
    for (unsigned i = 0; i < k; i++) {
        counter.count((kmer >> ((k-i-1)*2)) & 3);
    }
}

void skip(SpecificCompositionCounter &counter, unsigned window) {
    for (unsigned i = 0; i < window; i++) {
        counter.count(0);
    }
}

void test_counter(SpecificCompositionCounter &counter,
                  unsigned k, unsigned window, unsigned uPatternKmer, unsigned uKmer) {
    // 0, RC(0) - kmers that mustn't be counted
    // 1, RC(1) - kmers that must be counted

    // Pattern at borders and RC in the middle, no overlap
    counter.initGene(0);

    unsigned flank = 20;
    put_kmer(counter, k, uPatternKmer);
    put_kmer(counter, k, Utils::reverseComplement(uKmer, k));
    skip(counter, flank);
    put_kmer(counter, k, uKmer);
    put_kmer(counter, k, Utils::reverseComplement(uPatternKmer, k));
    put_kmer(counter, k, Utils::reverseComplement(uKmer, k));
    skip(counter, flank);
    put_kmer(counter, k, uKmer);
    put_kmer(counter, k, uPatternKmer);

    counter.finalizeGene();

    // Empty gene
    counter.initGene(1);
    counter.finalizeGene();

    // border patterns, windows overlap with patterns
    counter.initGene(2);

    put_kmer(counter, k, uPatternKmer);
    skip(counter, window-1);
    put_kmer(counter, k, Utils::reverseComplement(uPatternKmer, k));

    counter.finalizeGene();
}

void test_counter2(SpecificCompositionCounter &counter,
                  unsigned k, unsigned window, unsigned uPatternKmer,
                  unsigned uKmer) {
    // 0, RC(0) - kmers that mustn't be counted
    // 1, RC(1) - kmers that must be counted

    // Pattern at borders and RC in the middle, no overlap
    counter.initGene(0);

    unsigned flank = 20;
    put_kmer(counter, k, uPatternKmer);
    put_kmer(counter, k, Utils::reverseComplement(uKmer, k));
    skip(counter, flank);
    put_kmer(counter, k, uKmer);
    put_kmer(counter, k, Utils::reverseComplement(uPatternKmer, k));
    put_kmer(counter, k, Utils::reverseComplement(uKmer, k));
    skip(counter, flank);
    put_kmer(counter, k, uKmer);
    put_kmer(counter, k, uPatternKmer);

    counter.finalizeGene();

    // Empty gene
    counter.initGene(1);
    counter.finalizeGene();

    // border patterns, windows overlap with patterns
    counter.initGene(2);
    put_kmer(counter, k, uPatternKmer);
    skip(counter, window-1);
    put_kmer(counter, k, Utils::reverseComplement(uPatternKmer, k));

    counter.finalizeGene();
}

context("SpecificCompositionCounter") {
    test_that("initialize") {
        std::string patternKmer("GCCG");
        std::string kmer("GGGG");
        unsigned window = 4;
        unsigned k = 4;
        unsigned uPatternKmer = Utils::stringToInt(patternKmer);
        unsigned uPatternKmerRC = Utils::reverseComplement(uPatternKmer, k);
        unsigned uKmer = Utils::stringToInt(kmer);
        unsigned uKmerRC = Utils::reverseComplement(uKmer, k);

        DataStructureFactory factory;

        std::vector<std::string> geneNames({"gene1", "gene2", "gene3"});

        auto generate_match = [k, window](unsigned spacer, unsigned id){
            return [k, window, spacer, id](unsigned a, unsigned b){
                unsigned mask = ((1 << (2*k)) - 1);
                unsigned masked_a = a & mask;
                bool result = ((a >> 2*k) == spacer + id*window) && (masked_a != 0) && (masked_a != mask);
                // if (!result) {
                //     std::cout << "MISMATCH!\n expected spacer:" << spacer
                //               << "\n expected id:" << id
                //               << "\n real spacer:" << (a >> 2*k) % window
                //               << "\n real id:" << (a >> 2*k) / window
                //               << "\n";
                // }
                return result;
            };
        };

        test_that("counter works with rc flag") {
            SpecificCompositionCounter counter(
                    factory, geneNames, std::vector<std::string>({patternKmer}),
                    k, 0, window-1,  true
            );
            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            Fake(Method(spy, sElementInput));

            test_counter(counter, k, window, uPatternKmer, uKmer);

            unsigned hexamer_count = Utils::getKmersTotal(k);
            CATCH_CHECK_NOTHROW(Verify(
                    Method(spy, sGeneInput).Using(0) +
                        Method(spy, sElementInput).Matching(generate_match(0, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(3, 0)) +

                        Method(spy, sElementInput).Matching(generate_match(3, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(0, 0)) +

                        Method(spy, sElementInput).Matching(generate_match(0, 1)) +
                        Method(spy, sElementInput).Matching(generate_match(1, 1)) +
                        Method(spy, sElementInput).Matching(generate_match(2, 1)) +
                        Method(spy, sElementInput).Matching(generate_match(3, 1)) +

                        Method(spy, sElementInput).Matching(generate_match(3, 1)) +
                        Method(spy, sElementInput).Matching(generate_match(2, 1)) +
                        Method(spy, sElementInput).Matching(generate_match(1, 1)) +
                        Method(spy, sElementInput).Matching(generate_match(0, 1)) +

                        Method(spy, sGeneInput).Using(1) +

                        Method(spy, sGeneInput).Using(2) +
                        Method(spy, sElementInput).Matching(generate_match(0, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(3, 0)) +

                        //Method(spy, sElementInput).Matching(generate_match(3, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match(0, 0))
            ));

            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(Method(spy, sElementInput)));
        }

        test_that("counter works with rc and multiple patterns flag") {
            //std::cout << "Multiple patterns\n";

            patternKmer = "TGTCTG";
            std::string patternKmer2("TGTCNN");
            kmer = "GGGGGG";
            k = 6;
            uPatternKmer = Utils::stringToInt(patternKmer);
            uPatternKmerRC = Utils::reverseComplement(uPatternKmer, k);
            uKmer = Utils::stringToInt(kmer);
            uKmerRC = Utils::reverseComplement(uKmer, k);

            auto generate_match2 = [k, window](unsigned spacer, unsigned id){
                return [k, window, spacer, id](unsigned a, unsigned b){
                    unsigned mask = ((1 << (2*k)) - 1);
                    unsigned masked_a = a & mask;
                    //std::cout << "Called: " << a << ' ' << b << '\n';
                    //std::cout << "got: " << (a >> 2*k) << ", expected: " << spacer + id*window << "\n";
                    return ((a >> 2*k) == spacer + id*window) && (masked_a != 0) && (masked_a != mask);
                };
            };

            SpecificCompositionCounter counter(
                    factory, geneNames,
                    std::vector<std::string>({patternKmer2, patternKmer}),
                    k, 0, window-1,  true
            );
            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            When(Method(spy, sElementInput)).AlwaysDo(
                [](unsigned a, unsigned b){}//std::cout << "sElementInput(" << a << ", " << b << ")\n";}
            );

            test_counter(counter, k, window, uPatternKmer, uKmer);

            unsigned hexamer_count = Utils::getKmersTotal(k);

            CATCH_CHECK_NOTHROW(Verify(
                    Method(spy, sGeneInput).Using(0) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(3, 0)) +

                        Method(spy, sElementInput).Matching(generate_match2(0, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(3, 1)) +

                        Method(spy, sElementInput).Matching(generate_match2(3, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 0)) +

                        Method(spy, sElementInput).Matching(generate_match2(0, 2)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 2)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 2)) +
                        Method(spy, sElementInput).Matching(generate_match2(3, 2)) +

                        Method(spy, sElementInput).Matching(generate_match2(3, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 1)) +

                        Method(spy, sElementInput).Matching(generate_match2(0, 3)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 3)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 3)) +
                        Method(spy, sElementInput).Matching(generate_match2(3, 3)) +

                        Method(spy, sElementInput).Matching(generate_match2(3, 2)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 2)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 2)) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 2)) +

                        Method(spy, sElementInput).Matching(generate_match2(3, 3)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 3)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 3)) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 3)) +

                        Method(spy, sGeneInput).Using(1) +

                        Method(spy, sGeneInput).Using(2) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(3, 0)) +

                        Method(spy, sElementInput).Matching(generate_match2(0, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(3, 1)) +

                        //Method(spy, sElementInput).Matching(generate_match2(3, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 0)) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 0)) +

                        //Method(spy, sElementInput).Matching(generate_match2(3, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(2, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(1, 1)) +
                        Method(spy, sElementInput).Matching(generate_match2(0, 1))
            ));
            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(Method(spy, sElementInput)));
        }

        test_that("counter works with no flags") {
            SpecificCompositionCounter counter(
                    factory, geneNames, std::vector<std::string>({patternKmer}),
                    k, 0, window-1, false
            );

            Mock<IDataStructure> spy(*(counter.getResult()));
            Fake(Method(spy, sGeneInput));
            Fake(Method(spy, sElementInput));

            test_counter(counter, k, window, uPatternKmer, uKmer);

            unsigned hexamer_count = Utils::getKmersTotal(k);
            CATCH_CHECK_NOTHROW(Verify(
                Method(spy, sGeneInput).Using(0) +
                Method(spy, sElementInput).Matching(generate_match(0, 0)) +
                Method(spy, sElementInput).Matching(generate_match(1, 0)) +
                Method(spy, sElementInput).Matching(generate_match(2, 0)) +
                Method(spy, sElementInput).Matching(generate_match(3, 0)) +

                Method(spy, sElementInput).Matching(generate_match(3, 1)) +
                Method(spy, sElementInput).Matching(generate_match(2, 1)) +
                Method(spy, sElementInput).Matching(generate_match(1, 1)) +
                Method(spy, sElementInput).Matching(generate_match(0, 1)) +

                Method(spy, sGeneInput).Using(1) +

                Method(spy, sGeneInput).Using(2) +
                Method(spy, sElementInput).Matching(generate_match(0, 0)) +
                Method(spy, sElementInput).Matching(generate_match(1, 0)) +
                Method(spy, sElementInput).Matching(generate_match(2, 0)) +
                Method(spy, sElementInput).Matching(generate_match(3, 0))
            ));
            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(Method(spy, sElementInput)));
        }
    }

    test_that("initialize fuzzy") {
        std::string patternKmer("GCCG");
        std::string kmer("GGGG");
        unsigned window = 4;
        unsigned k = 4;
        unsigned uPatternKmer = Utils::stringToInt(patternKmer);
        unsigned uKmer = Utils::stringToInt(kmer);
        unsigned uKmerRC = Utils::reverseComplement(uKmer, k);
        unsigned uPatternKmerRC = Utils::reverseComplement(uPatternKmer, k);
        unsigned mask = (1 << (2*k)) - 1;
        std::unique_ptr<Pattern> pattern(new Pattern(patternKmer));
        std::vector<std::string> geneNames({"gene1", "gene2", "gene3"});

        DataStructureFactory factory;

        test_that("counter works with fuzzy flags") {
            SpecificCompositionCounter counter(
                    factory, geneNames, std::vector<std::string>({patternKmer}),
                    k, 0, window-1, false, true, true, true
            );
            Mock<IDataStructure> data(*(counter.getResult()));
            Fake(Method(data, sGeneInput));
            Fake(Method(data, sElementInput));

            test_counter2(counter, k, window, uPatternKmer, uKmer);
            CATCH_CHECK_NOTHROW(Verify(Method(data, sGeneInput).Using(0) +
                Method(data, sElementInput).Using(uKmerRC, _) +
                Method(data, sElementInput).Using((uKmerRC << 2) & mask, _) +
                Method(data, sElementInput).Using((uKmerRC << 4) & mask, _) +
                Method(data, sElementInput).Using((uKmerRC << 6) & mask, _) +
                Method(data, sElementInput).Using(uKmer >> 6, _) +
                Method(data, sElementInput).Using(uKmer >> 4, _) +
                Method(data, sElementInput).Using(uKmer >> 2, _) +
                Method(data, sElementInput).Using(uKmer, _) +
                Method(data, sGeneInput).Using(1) +
                Method(data, sGeneInput).Using(2) +
                Method(data, sElementInput).Using(uPatternKmerRC >> 6, _) +
                Method(data, sElementInput).Using(uPatternKmerRC >> 4, _) +
                Method(data, sElementInput).Using(uPatternKmerRC >> 2, _) +
                Method(data, sElementInput).Using(uPatternKmerRC, _)));
            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(data));
        }

        test_that("counter works with rc flag") {
            pattern -> add("CGGC");
            SpecificCompositionCounter counter(
                    factory, geneNames, std::vector<std::string>({patternKmer}),
                    k, 0, window-1, true, true, true, true
            );

            Mock<IDataStructure> data(*(counter.getResult()));
            Fake(Method(data, sGeneInput));
            Fake(Method(data, sElementInput));

            test_counter2(counter, k, window, uPatternKmer, uKmer);

            auto matcher = [mask](unsigned a, unsigned b){
                bool result = (a != 0) && (a != mask);
                if (!result) {
                    // std::cout << "Matcher failure: " << a << ' ' << b;
                }
                return result;
            };
            CATCH_CHECK_NOTHROW(Verify(Method(data, sGeneInput).Using(0) +
                Method(data, sElementInput).Matching(matcher) * window * 4 +
                Method(data, sGeneInput).Using(1) +
                Method(data, sGeneInput).Using(2) +
                Method(data, sElementInput).Matching(matcher) * (window * 2-1)
            ));
            CATCH_CHECK_NOTHROW(VerifyNoOtherInvocations(Method(data, sElementInput)));
        }
    }
}
