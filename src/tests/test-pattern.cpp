#include <testthat.h>
#include <iostream>

#include "../Pattern/Pattern.h"
#include "../Pattern/PatternWrongSizeException.hpp"
#include "../Utils/Utils.h"

#define K_ 5
#define SIZE (1 << K_ * 2)

void testAllKmers(bool * testData, Pattern& pattern, unsigned size) {
    for (unsigned i = 0; i < size; i++) {
        expect_true(pattern.check(i) == testData[i]);
    }
}

void markKmers(std::string * kmers, unsigned kmersCount,
               bool * testData, unsigned size) {
    for (unsigned i = 0; i< size; i++) {
        testData[i] = false;
    }
    for (unsigned i = 0; i< kmersCount; i++) {
        testData[Utils::stringToInt(kmers[i])] = true;
    }
}

void test_find(const std::string& kmerString, const std::vector<unsigned> & kmers) {
    unsigned kmer = Utils::stringToInt(kmerString);
    expect_true(*std::find(kmers.begin(), kmers.end(), kmer) == kmer);
}

context ("Pattern") {
    test_that("string constructor works") {
        std::string p = std::string("ACACA");
        Pattern pattern(p);

        bool testData[SIZE];

        markKmers(&p, 1, testData, SIZE);
        testAllKmers(testData, pattern, SIZE);
    }

    test_that("empty constructor + add works") {
        std::string p = std::string("ACACA");
        Pattern pattern;
        pattern.add(p);

        bool testData[SIZE];

        markKmers(&p, 1, testData, SIZE);
        testAllKmers(testData, pattern, SIZE);
    }

    test_that("ambigious pattern translates into all possible unambigious ones") {
        std::string p = std::string("ACANB");
        Pattern pattern(p);
        std::string kmers[12] = {
            "ACAAC", "ACAAG", "ACAAT",
            "ACACC", "ACACG", "ACACT",
            "ACAGC", "ACAGG", "ACAGT",
            "ACATC", "ACATG", "ACATT",
        };
        bool testData[SIZE];
        markKmers(kmers, 12, testData, SIZE);
        testAllKmers(testData, pattern, SIZE);
    }

    test_that("adding longer kmer throws exception") {
        std::string p = std::string("ACANB");
        Pattern pattern(p);
        expect_error_as(pattern.add("TTTTTT"), PatternWrongSizeException);
    }

    test_that("adding shorter kmer throws exception") {
        std::string p = std::string("ACANB");
        Pattern pattern(p);
        expect_error_as(pattern.add("TTTT"), PatternWrongSizeException);
    }

    test_that("getKmers works") {
        std::string p = std::string("ACANB");
        Pattern pattern(p);
        std::vector<unsigned> kmers = pattern.getKmers();
        expect_true(kmers.size() == 12);
        std::string answers[12] = {
            "ACAAC", "ACAAG", "ACAAT",
            "ACACC", "ACACG", "ACACT",
            "ACAGC", "ACAGG", "ACAGT",
            "ACATC", "ACATG", "ACATT",
        };
        for (std::string kmer : answers) {
            test_find(kmer, kmers);
        }
    }

    test_that("copy constructor works")
    {
        // Copy of ambigious_pattern
        std::string p = std::string("ACANB");
        Pattern pattern(p);
        Pattern pattern2(pattern);
        std::string kmers[12] = {
            "ACAAC", "ACAAG", "ACAAT",
            "ACACC", "ACACG", "ACACT",
            "ACAGC", "ACAGG", "ACAGT",
            "ACATC", "ACATG", "ACATT",
        };
        bool testData[SIZE];
        markKmers(kmers, 12, testData, SIZE);
        testAllKmers(testData, pattern2, SIZE);
    }
}
