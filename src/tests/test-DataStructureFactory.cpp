#include <testthat.h>
#include <iostream>

#include "../DataStructures/DataStructureFactory.h"
#include "../DataStructures/ElementCounts.h"
#include "../DataStructures/MotifPositionsSparse.h"

context("DataStructureFactory") {
    test_that("initialization") {
        DataStructureFactory factory;
        factory.setType(DataStructureFactory::type::ElementCounts);
        expect_true(factory.getType() == DataStructureFactory::type::ElementCounts);

        test_that("copy works") {
            DataStructureFactory another(factory);
            expect_true(another.getType() == DataStructureFactory::type::ElementCounts);
        }
    }

    test_that("structures are constructed with correct types") {
        std::vector<std::string> elementLabels({"elem1", "elem2", "elem3"});
        std::vector<std::string> geneLabels({"gene1", "gene2", "gene3", "gene4"});

        DataStructureFactory factory;
        factory.setType(DataStructureFactory::type::ElementCounts);

        std::function<std::string(unsigned)> labelGenerator =
            [](unsigned id){return "elem" + std::to_string(id);};
        IDataStructure * data = factory.create(labelGenerator, geneLabels);
        expect_true(dynamic_cast<ElementCounts *>(data) != 0);
        delete data;

        factory.setType(DataStructureFactory::type::MotifPositionsSparse);
        data = factory.create(labelGenerator, geneLabels);
        expect_true(dynamic_cast<MotifPositionsSparse *>(data) != 0);
        delete data;
    }
}
