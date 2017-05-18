#include <testthat.h>
#include <iostream>

#include "../Motifs/encodings.h"

context("encodings") {
    test_that ("put IUPAC works") {
        cell code = 0;
        putIUPAC(&code, 0, 7);
        expect_true(code == 7);

        putIUPAC(&code, 0, 8);
        expect_true(code == 8);

        putIUPAC(&code, 2, 8);
        expect_true(code == 0x808);

        cell longcode[2];
        longcode[0] = longcode[1] = 0;
        putIUPAC(longcode, 18, 0xf);

        expect_true(longcode[0] == 0);
        expect_true(longcode[1] == 0xf00);
    }

    test_that ("get IUPAC works") {
        cell longcode[2];
        longcode[0] = 0x1248;
        longcode[1] = 0x1248;

        expect_true(getIUPAC(longcode, 0) == 8);
        expect_true(getIUPAC(longcode, 1) == 4);
        expect_true(getIUPAC(longcode, 2) == 2);
        expect_true(getIUPAC(longcode, 3) == 1);

        expect_true(getIUPAC(longcode, 16) == 8);
        expect_true(getIUPAC(longcode, 17) == 4);
        expect_true(getIUPAC(longcode, 18) == 2);
        expect_true(getIUPAC(longcode, 19) == 1);
    }

    test_that ("put Compact works") {
        cell code = 0;
        putCompact(&code, 0, 3);
        expect_true(code == 3);

        putCompact(&code, 0, 1);
        expect_true(code == 1);

        putCompact(&code, 2, 2);
        expect_true(code == 0x21);

        cell longcode[2];
        longcode[0] = longcode[1] = 0;
        putCompact(longcode, 34, 3);

        expect_true(longcode[0] == 0);
        expect_true(longcode[1] == 0x30);
    }

    test_that ("get compact works") {
        cell longcode[2];
        longcode[0] = 0xe4;
        longcode[1] = 0xe4;

        expect_true(getCompact(longcode, 0) == 0);
        expect_true(getCompact(longcode, 1) == 1);
        expect_true(getCompact(longcode, 2) == 2);
        expect_true(getCompact(longcode, 3) == 3);

        expect_true(getCompact(longcode, 32) == 0);
        expect_true(getCompact(longcode, 33) == 1);
        expect_true(getCompact(longcode, 34) == 2);
        expect_true(getCompact(longcode, 35) == 3);
    }
}
