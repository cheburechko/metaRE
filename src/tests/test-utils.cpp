#include <testthat.h>
#include <iostream>

#include "../Utils/Utils.h"

using namespace std;
context("Utils works") {
    test_that ("reverseComplementChar works") {
        expect_true(Utils::reverseComplementChar('a') == 't');
        expect_true(Utils::reverseComplementChar('t') == 'a');
        expect_true(Utils::reverseComplementChar('c') == 'g');
        expect_true(Utils::reverseComplementChar('g') == 'c');

        expect_true(Utils::reverseComplementChar('A', true) == 'T');
        expect_true(Utils::reverseComplementChar('T', true) == 'A');
        expect_true(Utils::reverseComplementChar('G', true) == 'C');
        expect_true(Utils::reverseComplementChar('C', true) == 'G');

        expect_error_as(Utils::reverseComplementChar('x'), UtilsNucleotideException);
    }

    test_that ("reverseComplementInt works") {
        expect_true(Utils::reverseComplementInt(0) == 3);
        expect_true(Utils::reverseComplementInt(1) == 2);
        expect_true(Utils::reverseComplementInt(2) == 1);
        expect_true(Utils::reverseComplementInt(3) == 0);

        expect_error_as(Utils::reverseComplementInt(10), UtilsNucleotideException);
    }

    test_that ("intToChar works") {
        expect_true(Utils::intToChar(0) == 'a');
        expect_true(Utils::intToChar(1) == 'c');
        expect_true(Utils::intToChar(2) == 'g');
        expect_true(Utils::intToChar(3) == 't');

        expect_true(Utils::intToChar(0, true) == 'A');
        expect_true(Utils::intToChar(1, true) == 'C');
        expect_true(Utils::intToChar(2, true) == 'G');
        expect_true(Utils::intToChar(3, true) == 'T');

        expect_error_as(Utils::intToChar(4), UtilsNucleotideException);
    }

    test_that("charToInt works") {
        expect_true(Utils::charToInt('a') == 0);
        expect_true(Utils::charToInt('c') == 1);
        expect_true(Utils::charToInt('g') == 2);
        expect_true(Utils::charToInt('t') == 3);

        expect_true(Utils::charToInt('A') == 0);
        expect_true(Utils::charToInt('C') == 1);
        expect_true(Utils::charToInt('G') == 2);
        expect_true(Utils::charToInt('T') == 3);

        expect_true(Utils::charToInt('x') == -1);
    }

    test_that("reverseComplementIntString works") {
        // a -> t
        expect_true(Utils::reverseComplement(0, 1) == 3);
        // ac -> gt
        expect_true(Utils::reverseComplement(1, 2) == 11);
        // taag -> ctta
        expect_true(Utils::reverseComplement(194, 4) == 124);

        expect_error_as(Utils::reverseComplement(350, 4), UtilsKmerTooBigException);
    }

    test_that("reverseComplementCharString works") {
        expect_true(Utils::reverseComplement(string("a")).compare("t") == 0);
        expect_true(Utils::reverseComplement(string("ac")).compare("gt") == 0);
        expect_true(Utils::reverseComplement(string("taag")).compare("ctta") == 0);
        expect_true(Utils::reverseComplement(string("taag"), true).compare("CTTA") == 0);
        expect_true(Utils::reverseComplement(string("TAAG")).compare("ctta") == 0);
    }

    test_that("stringToInt works") {
        expect_true(Utils::stringToInt("taag") == 194);
        expect_true(Utils::stringToInt("actg") == 30);
    }

    test_that("intToString works") {
        expect_true(Utils::intToString(194, 4).compare("taag") == 0);
        expect_true(Utils::intToString(30, 4).compare("actg") == 0);
        expect_true(Utils::intToString(30, 4, true).compare("ACTG") == 0);

        expect_error_as(Utils::intToString(350, 4), UtilsKmerTooBigException);
    }

    test_that("getKmersTotal works") {
        expect_true(Utils::getKmersTotal(6) == 4096);
    }

    test_that("allReverseComplements works") {
        unsigned k = 3;
        auto rcs = Utils::allReverseComplements(k);
        // cgg -> ccg
        expect_true(rcs[26] == 22);
        // aaa -> ttt
        expect_true(rcs[0] == 63);
        // att -> aat
        expect_true(rcs[15] == 3);
    }

    test_that("getPlusMinusMapping works") {
        unsigned k = 3;
        auto pm = Utils::getPlusMinusMapping(k);
        // cgg -> ccg
        expect_true(pm[26] == pm[22]);
        // att -> aat
        expect_true(pm[15] == pm[3]);
    }

    test_that("generateAllKmerLabels works") {
        unsigned k = 4;
        std::vector<std::string> labels = Utils::generateAllKmerLabels(k);
        for (unsigned i = 0; i < Utils::getKmersTotal(k); i++) {
            std::string kmer = Utils::intToString(i, k);
            expect_true(std::find(labels.begin(), labels.end(), kmer) != labels.end());
        }
    }

}
