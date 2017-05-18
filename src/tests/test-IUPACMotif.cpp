#include <testthat.h>
#include <iostream>

#include "../Motifs/IUPACMotif.h"

context("IUPACMotif") {
    test_that("initialize") {
        IUPACMotif motif("ACGT");
        IUPACMotif longMotif("AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT");
        IUPACMotif degenMotif("BAD-WYSDM");
        IUPACMotif longDegenMotif("AAAACCCCNGGGTTTTAAAACCCCGGGGTTTT");

        test_that("degenerate flag is correct") {
            expect_false(motif.isDegenerate());
            expect_true(degenMotif.isDegenerate());
            expect_false(longMotif.isDegenerate());
            expect_true(longDegenMotif.isDegenerate());
        }

        test_that("length is correct") {
            expect_true(motif.getLength() == 4);
            expect_true(degenMotif.getLength() == 9);
            expect_true(IUPACMotif("AAAAAAAAAAAA").getLength() == 12);
            expect_true(IUPACMotif("AAAAAAAAAAAAAAAAAAAAAAAA").getLength() == 24);
            expect_true(IUPACMotif("AAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCC").getLength() == 48);
        }

        test_that("motif encodes to IUPAC with array generation") {
            cell * code = motif.encodeIUPAC();
            expect_true(*code == 0x1248);
            delete code;

            code = degenMotif.encodeIUPAC();
            expect_true(*code == 0xe1d09a6d3);
            delete code;

            test_that("long motifs also encode") {
                cell * code = longMotif.encodeIUPAC();
                expect_true(code[3] == 0x1111111111111111);
                expect_true(code[2] == 0x2222222222222222);
                expect_true(code[1] == 0x4444444444444444);
                expect_true(code[0] == 0x8888888888888888);
                delete code;

                code = longDegenMotif.encodeIUPAC();
                expect_true(code[1] == 0x11112222f4448888);
                expect_true(code[0] == 0x1111222244448888);
                delete code;
            }
        }

        test_that("motif encodes to IUPAC in preset array") {
            cell code = 0;
            motif.encodeIUPAC(&code);
            expect_true(code == 0x1248);

            code = 0;
            degenMotif.encodeIUPAC(&code);
            expect_true(code == 0xe1d09a6d3);

            test_that("long motifs also encode") {
                cell longCode[4];
                longMotif.encodeIUPAC(longCode);
                expect_true(longCode[3] == 0x1111111111111111);
                expect_true(longCode[2] == 0x2222222222222222);
                expect_true(longCode[1] == 0x4444444444444444);
                expect_true(longCode[0] == 0x8888888888888888);

                longDegenMotif.encodeIUPAC(longCode);
                expect_true(longCode[1] == 0x11112222f4448888);
                expect_true(longCode[0] == 0x1111222244448888);
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

        test_that("degenrate motifs can't be encoded in compact") {
            expect_error(degenMotif.encodeCompact());
            expect_error(longDegenMotif.encodeCompact());
        }

        test_that("motif encodes to compact form with array generation") {
            cell code = 0;
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
            expect_true(degenMotif.getString().compare("bad-wysdm") == 0);
            expect_true(longDegenMotif.getString().compare("aaaaccccngggttttaaaaccccggggtttt") == 0);
        }
    };

    test_that("initialize RC/include cases") {
        IUPACMotif motif("AAGG");
        IUPACMotif longMotif("AAACCTTGGAATTTGCGC");
        IUPACMotif degenMotif("-ACGTRYKMSWBDHVN");

        test_that("reverse complements are generated correctly") {
            Motif * bufm = motif.getReverseComplement();
            std::string buf = bufm -> getString();
            expect_true(buf.compare("cctt") == 0);
            delete bufm;

            bufm = longMotif.getReverseComplement();
            buf = bufm -> getString();
            expect_true(buf.compare("gcgcaaattccaaggttt") == 0);
            delete bufm;

            bufm = degenMotif.getReverseComplement();
            buf = bufm -> getString();
            expect_true(buf.compare("nbdhvwskmryacgt-") == 0);
            delete bufm;
        }

        test_that("motif recognizes its reverse complement") {
            expect_true(motif.isReverseComplementOf(IUPACMotif("CCTT")));
            expect_false(motif.isReverseComplementOf(IUPACMotif("CCTG")));

            IUPACMotif buf("GCGCAAATTCCAAGGTTT");
            Motif * bufRC = buf.getReverseComplement();
            expect_true( longMotif.isReverseComplementOf(buf));
            expect_false(longMotif.isReverseComplementOf(IUPACMotif("GCGCAAATTGCAAGGTTT")));
            delete bufRC;

            IUPACMotif buf2("nbdhvwskmryacgt-");
            bufRC = buf2.getReverseComplement();
            expect_true( degenMotif.isReverseComplementOf(buf2));
            expect_false(degenMotif.isReverseComplementOf(IUPACMotif("nbdhvwstmryacgt-")));
            delete bufRC;
        }

        test_that("inclusion is checked correctly") {
            expect_true(motif.includes(motif));
            expect_false(motif.includes(IUPACMotif("AAAC")));
            expect_false(motif.includes(IUPACMotif("TGTC")));
            expect_false(motif.includes(IUPACMotif("AAG")));
            expect_false(motif.includes(IUPACMotif("AAAGG")));

            expect_true(longMotif.includes(longMotif));
            expect_false(longMotif.includes(IUPACMotif("AAACCGTGGAATTTGCGC")));
            expect_false(longMotif.includes(IUPACMotif("TGTGTCGAGCAGTGCATG")));
            expect_false(longMotif.includes(IUPACMotif("AAACCGTGGAATTTGCG")));
            expect_false(longMotif.includes(IUPACMotif("AAACCGTGGAATTTGCGCC")));

            IUPACMotif degenMotif("ACGTRYKMSWBDHVN");
            expect_true(degenMotif.includes(degenMotif));
            expect_true(degenMotif.includes(IUPACMotif("acgtacgacacaccc")));
            expect_false(degenMotif.includes(IUPACMotif("aacgtacgacacaccc")));
            expect_false(degenMotif.includes(IUPACMotif("aggtacgacacaccc")));
            expect_false(degenMotif.includes(IUPACMotif("acgtacgacaaaccc")));
            expect_false(degenMotif.includes(IUPACMotif("acgttcgacacaccc")));
        }

        test_that("equivalence is checked correctly") {
            expect_true(motif == IUPACMotif("AAGG"));
            expect_false(motif == IUPACMotif("AAAC"));
            expect_false(motif == IUPACMotif("TGTC"));
            expect_false(motif == IUPACMotif("AAG"));
            expect_false(motif == IUPACMotif("AAAGG"));

            expect_true(longMotif == IUPACMotif("AAACCTTGGAATTTGCGC"));
            expect_false(longMotif == IUPACMotif("AAACCGTGGAATTTGCGC"));
            expect_false(longMotif == IUPACMotif("TGTGTCGAGCAGTGCATG"));
            expect_false(longMotif == IUPACMotif("AAACCGTGGAATTTGCG"));
            expect_false(longMotif == IUPACMotif("AAACCGTGGAATTTGCGCC"));

            expect_true(degenMotif == IUPACMotif("-ACGTRYKMSWBDHVN"));
            expect_false(degenMotif == IUPACMotif("--ACGTRYKMSWBDHVN"));
            expect_false(degenMotif == IUPACMotif("-ACGTRYKMSWBDHVN-"));
            expect_false(degenMotif == IUPACMotif("ACGTRYKMSWBDHVN"));
            expect_false(degenMotif == IUPACMotif("-ACGTRYKMSWBDHV"));
            expect_false(degenMotif == IUPACMotif("-ACGTRYTMSWBDHVN"));

        }
    };
};
