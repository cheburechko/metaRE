#include <testthat.h>
#include <iostream>

#include "../Motifs/CompactMotif.h"

context("CompactMotif") {
    test_that("initialize") {
        CompactMotif motif("ACGT");
        CompactMotif longMotif("AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT");

        test_that("length is correct") {
            expect_true(motif.getLength() == 4);
            expect_true(CompactMotif("AAAAAAAAAAAA").getLength() == 12);
            expect_true(CompactMotif("AAAAAAAAAAAAAAAAAAAAAAAA").getLength() == 24);
            expect_true(CompactMotif("AAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCC").getLength() == 48);
        }

        test_that("motif encodes to IUPAC with array generation") {
            cell * code = motif.encodeIUPAC();
            expect_true(*code == 0x1248);
            delete code;

            test_that("long motifs also encode") {
                cell * code = longMotif.encodeIUPAC();
                expect_true(code[3] == 0x1111111111111111);
                expect_true(code[2] == 0x2222222222222222);
                expect_true(code[1] == 0x4444444444444444);
                expect_true(code[0] == 0x8888888888888888);
                delete code;
            }
        }

        test_that("motif encodes to IUPAC in preset array") {
            cell code = 0;
            motif.encodeIUPAC(&code);
            expect_true(code == 0x1248);

            test_that("long motifs also encode") {
                cell longCode[4];
                longMotif.encodeIUPAC(longCode);
                expect_true(longCode[3] == 0x1111111111111111);
                expect_true(longCode[2] == 0x2222222222222222);
                expect_true(longCode[1] == 0x4444444444444444);
                expect_true(longCode[0] == 0x8888888888888888);
            }
        }

        test_that("motif encodes to compact form with array generation") {
            cell * code = motif.encodeCompact();
            expect_true(*code == 0b00011011);
            delete code;

            test_that("long motif also encode") {
                code = longMotif.encodeCompact();
                expect_true(code[1] == 0x0000000055555555);
                expect_true(code[0] == 0xaaaaaaaaffffffff);
                delete code;
            }
        }

        test_that("motif encodes to compact form with array generation") {
            cell code;
            motif.encodeCompact(&code);
            expect_true(code == 0b00011011);

            test_that("long motif also encode") {
                cell longCode[2];
                longMotif.encodeCompact(longCode);
                expect_true(longCode[1] == 0x0000000055555555);
                expect_true(longCode[0] == 0xaaaaaaaaffffffff);
            }
        }

        test_that("string is generated") {
            expect_true(motif.getString().compare("acgt") == 0);
            expect_true(longMotif.getString().compare("aaaaaaaaaaaaaaaaccccccccccccccccggggggggggggggggtttttttttttttttt") == 0);
        }

    };

    test_that("initialize RC/include cases") {
        CompactMotif motif("AAGG");
        CompactMotif longMotif("AAACCTTGGAATTTGCGC");

        test_that("reverse complements are generated correctly") {
            Motif * bufm = motif.getReverseComplement();
            std::string buf = bufm -> getString();
            expect_true(buf.compare("cctt") == 0);
            delete bufm;

            bufm = longMotif.getReverseComplement();
            buf = bufm -> getString();
            expect_true(buf.compare("gcgcaaattccaaggttt") == 0);
            delete bufm;
        }

        test_that("motif recognizes its reverse complement") {
            expect_true(motif.isReverseComplementOf(CompactMotif("CCTT")));
            expect_false(motif.isReverseComplementOf(CompactMotif("CCTG")));

            CompactMotif buf("GCGCAAATTCCAAGGTTT");
            Motif * bufRC = buf.getReverseComplement();
            expect_true( longMotif.isReverseComplementOf(buf));
            expect_false(longMotif.isReverseComplementOf(CompactMotif("GCGCAAATTGCAAGGTTT")));
            delete bufRC;
        }

        test_that("inclusion is checked correctly") {
            expect_true(motif.includes(motif));
            expect_false(motif.includes(CompactMotif("AAAC")));
            expect_false(motif.includes(CompactMotif("TGTC")));
            expect_false(motif.includes(CompactMotif("AAG")));
            expect_false(motif.includes(CompactMotif("AAAGG")));

            expect_true(longMotif.includes(longMotif));
            expect_false(longMotif.includes(CompactMotif("AAACCGTGGAATTTGCGC")));
            expect_false(longMotif.includes(CompactMotif("TGTGTCGAGCAGTGCATG")));
            expect_false(longMotif.includes(CompactMotif("AAACCGTGGAATTTGCG")));
            expect_false(longMotif.includes(CompactMotif("AAACCGTGGAATTTGCGCC")));
        }

        test_that("equivalence is checked correctly") {
            expect_true(motif == CompactMotif("AAGG"));
            expect_false(motif == CompactMotif("AAAC"));
            expect_false(motif == CompactMotif("TGTC"));
            expect_false(motif == CompactMotif("AAG"));
            expect_false(motif == CompactMotif("AAAGG"));

            expect_true(longMotif == CompactMotif("AAACCTTGGAATTTGCGC"));
            expect_false(longMotif == CompactMotif("AAACCGTGGAATTTGCGC"));
            expect_false(longMotif == CompactMotif("TGTGTCGAGCAGTGCATG"));
            expect_false(longMotif == CompactMotif("AAACCGTGGAATTTGCG"));
            expect_false(longMotif == CompactMotif("AAACCGTGGAATTTGCGCC"));
        }
    };
};
