#include <testthat.h>
#include <iostream>

#include "../DataStructures/ElementCounts.h"

context("ElementCounts") {
    test_that("initialization works") {
        std::vector<std::string> elementLabels({"elem0", "elem1"});
        std::vector<std::string> geneLabels({"gene1", "gene2", "gene3", "gene4"});
        std::function<std::string(unsigned)> labelGenerator =
            [](unsigned id){return "elem" + std::to_string(id);};
        ElementCounts data(labelGenerator, geneLabels);

        test_that("gene getters work") {
            expect_true(data.getGeneCount() == geneLabels.size());

            const std::vector<std::string> & dataGeneLabels = data.getGeneLabels();
            for (unsigned i = 0; i < geneLabels.size(); i++) {
                expect_true(data.getGeneLabel(i).compare(geneLabels[i]) == 0);
                expect_true(dataGeneLabels[i].compare(geneLabels[i]) == 0);
            }
        };
        test_that("structure getters work") {
            data.sGeneInput(0);
            data.sElementInput(0, 10);
            data.sElementInput(1, 20);
            data.sGeneInput(1);
            data.sElementInput(0, 11);
            data.sElementInput(0, 21);
            data.sElementInput(0, 31);
            data.sGeneInput(2);
            data.sElementInput(0, 12);
            data.sElementInput(1, 22);
            data.sGeneInput(0);
            data.sElementInput(0, 13);
            data.sElementInput(1, 23);

            const std::unordered_map<unsigned, std::vector<int> >& structure = data.getStructure();
            expect_true(structure.size() == 2);
            expect_true(structure.at(0).size() == 1);
            expect_true(structure.at(0)[0] == 6);
            expect_true(structure.at(1).size() == 1);
            expect_true(structure.at(1)[0] == 3);

            test_that("element getters work") {
                expect_true(data.getElementCount() == elementLabels.size());

                auto dataElementLabels = data.getElementLabels();
                for (unsigned i = 0; i < elementLabels.size(); i++) {
                    expect_true(data.getElementLabel(i).compare(elementLabels[i]) == 0);
                    expect_true(dataElementLabels.at(i).compare(elementLabels[i]) == 0);
                }
            };
        };
    };
}
