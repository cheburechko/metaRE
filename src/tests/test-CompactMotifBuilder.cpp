#include <testthat.h>
#include <iostream>

#include "../Motifs/CompactMotifBuilder.h"

context("CompactMotifBuilder") {
    test_that("initialize") {
        unsigned size = 4;
        CompactMotifBuilder builder(size);

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
            builder.put(CHAR_TO_COMPACT('a'));
            builder.put(CHAR_TO_COMPACT('c'));
            builder.put(CHAR_TO_COMPACT('g'));
            builder.put(CHAR_TO_COMPACT('t'));
            builder.put(CHAR_TO_COMPACT('g'));
            builder.put(CHAR_TO_COMPACT('c'));
            builder.put(CHAR_TO_COMPACT('t'));


            std::string answer("tgct");
            std::string complementAnswer("agca");
            for (unsigned i = 1; i <= size; i++ ) {
                expect_true(builder.ready(i));

                std::unique_ptr<CompactMotif> motif(builder.build(i));
                CATCH_INFO("Motif: " << motif->getString());
                CATCH_INFO("Substring: " << answer.substr(size-i, i));
                expect_true(answer.substr(size-i, i).compare(motif->getString()) == 0);

                std::unique_ptr<CompactMotif> complement(builder.buildComplement(i));
                CATCH_INFO("Complement: " << complement->getString());
                CATCH_INFO("Substring: " << complementAnswer.substr(0, i));
                expect_true(complementAnswer.substr(0, i).compare(complement->getString()) == 0);
            }

            cell buf = 0;
            builder.write(&buf, 3);
            expect_true(buf == 0x27);
        }

        test_that("skipping works") {
            for (unsigned i = 0; i < size; i++) {
                builder.put(0);
            }
            expect_true(builder.ready());
            builder.skip();
            std::string answer("cccccccccc");
            for (unsigned i = 1; i <= size; i++) {
                for (unsigned j = 1; j <= size; j++) {
                    expect_true(builder.ready(j) == (j < i));
                }
                builder.put(1);
                std::unique_ptr<CompactMotif> motif(builder.build(i));
                expect_true(answer.substr(size-i, i).compare(motif->getString()) == 0);
            }
        }

        test_that("put IUPAC works as expected") {
            builder.putIUPAC(1);
            builder.putIUPAC(2);
            builder.putIUPAC(4);
            builder.putIUPAC(8);

            expect_true(builder.ready(4));
            std::unique_ptr<CompactMotif> motif(builder.build(4));
            expect_true(motif->getString().compare("acgt")==0);

            expect_error(builder.putIUPAC(3));
        }

    }
}
