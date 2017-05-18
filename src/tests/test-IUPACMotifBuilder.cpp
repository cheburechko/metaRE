#include <testthat.h>
#include <iostream>

#include "../Motifs/IUPACMotifBuilder.h"

context("IUPACMotif") {
    test_that("initialize") {
        unsigned size = 10;
        IUPACMotifBuilder builder(size);

        test_that("builder is not ready until it's filled up") {
            for (unsigned input = 1; input <= 2*size; input++) {
                builder.put(1);
                for (unsigned test = 1; test <= size; test++) {
                    expect_true(builder.ready(test) == (test <= input));
                }
                expect_true(builder.ready() == (input >= size));
            }
        }

        test_that("builder creates motifs") {
            builder.put(CHAR_TO_IUPAC('a'));
            builder.put(CHAR_TO_IUPAC('c'));
            builder.put(CHAR_TO_IUPAC('g'));
            builder.put(CHAR_TO_IUPAC('t'));
            builder.put(CHAR_TO_IUPAC('r'));
            builder.put(CHAR_TO_IUPAC('y'));
            builder.put(CHAR_TO_IUPAC('m'));
            builder.put(CHAR_TO_IUPAC('k'));
            builder.put(CHAR_TO_IUPAC('n'));
            builder.put(CHAR_TO_IUPAC('d'));
            builder.put(CHAR_TO_IUPAC('b'));
            builder.put(CHAR_TO_IUPAC('v'));
            builder.put(CHAR_TO_IUPAC('-'));
            builder.put(CHAR_TO_IUPAC('h'));
            builder.put(CHAR_TO_IUPAC('s'));
            builder.put(CHAR_TO_IUPAC('w'));

            std::string answer("mkndbv-hsw");
            std::string complementAnswer("wsd-bvhnmk");
            for (unsigned i = 1; i <= size; i++ ) {
                expect_true(builder.ready(i));

                std::unique_ptr<IUPACMotif> motif(builder.build(i));
                CATCH_INFO("Motif: " << motif->getString());
                CATCH_INFO("Substring: " << answer.substr(size-i, i));
                expect_true(answer.substr(size-i, i).compare(motif->getString()) == 0);

                std::unique_ptr<IUPACMotif> complement(builder.buildComplement(i));
                CATCH_INFO("Complement: " << complement->getString());
                CATCH_INFO("Substring: " << complementAnswer.substr(0, i));
                expect_true(complementAnswer.substr(0, i).compare(complement->getString()) == 0);
            }
        }

        test_that("skipping works") {
            for (unsigned i = 0; i < size; i++) {
                builder.put(1);
            }
            expect_true(builder.ready());
            builder.skip();
            std::string answer("cccccccccc");
            for (unsigned i = 1; i <= size; i++) {
                for (unsigned j = 1; j <= size; j++) {
                    expect_true(builder.ready(j) == (j < i));
                }
                builder.put(2);
                std::unique_ptr<IUPACMotif> motif(builder.build(i));
                expect_true(answer.substr(size-i, i).compare(motif->getString()) == 0);
            }
        }

        test_that("put compact works as expected") {
            builder.putCompact(0);
            builder.putCompact(1);
            builder.putCompact(2);
            builder.putCompact(3);

            expect_true(builder.ready(4));
            std::unique_ptr<IUPACMotif> motif(builder.build(4));
            expect_true(motif->getString().compare("acgt")==0);
        }

        test_that("matching works") {
            builder.put(CHAR_TO_IUPAC('a'));
            builder.put(CHAR_TO_IUPAC('c'));
            builder.put(CHAR_TO_IUPAC('g'));
            builder.put(CHAR_TO_IUPAC('t'));
            builder.put(CHAR_TO_IUPAC('a'));
            builder.put(CHAR_TO_IUPAC('c'));
            builder.put(CHAR_TO_IUPAC('g'));
            builder.put(CHAR_TO_IUPAC('t'));

            IUPACMotif pattern("acgt");
            expect_true(builder.matches(pattern));
            IUPACMotif pattern2("NNSWDBRY");
            expect_true(builder.matches(pattern2));
            IUPACMotif pattern3("NNSWDBRV");
            expect_false(builder.matches(pattern3));
            IUPACMotif pattern4("acct");
            expect_false(builder.matches(pattern4));
            IUPACMotif pattern5("acgtac");
            expect_false(builder.matches(pattern5));
            expect_true(builder.matches(pattern5, true));
        }

        // TODO test the rest of the functions
    }
}
