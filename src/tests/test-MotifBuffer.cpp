#include <testthat.h>
#include <iostream>

#include "../Counters/MotifBuffer.hpp"

context("MotifBuffer") {
    test_that("initialize") {
        unsigned size = 11;
        MotifBuffer<cell> buffer(11);

        expect_true(buffer.size() == size);

        test_that("buffer is not ready until it's filled up") {
            int center = size / 2;
            expect_true(center == buffer.getCenterPos());
            for (unsigned input = 0; input < size; input++) {
                int value = buffer.put(input);
                expect_true(buffer.centerAvailable() == (input >= center));
            }
            expect_true(buffer.getCenter() == center);
            for (unsigned i = 0; i < size; i++) {
                expect_true(buffer.isValid(i));
                expect_true(buffer.get(i) == i);
                expect_true(buffer.absPosition(i) == i);
                expect_true(buffer.posFromCenter(i) ==  i-center);
            }

            test_that("buffer reset works") {
                buffer.reset();
                expect_false(buffer.centerAvailable());
                for (unsigned i = 0; i < size; i++) {
                    expect_false(buffer.isValid(i));
                }
            }

            test_that("putting pushes out the right values") {
                for (unsigned input = 0; input < size; input++) {
                    int value = buffer.put(size+1);
                    expect_true(input == value);
                }

                test_that("absPosition is computed correctly") {
                    for (unsigned input = 0; input < size; input++) {
                        expect_true(buffer.absPosition(input) == input+size);
                    }
                }
            }

            test_that("skipping works") {
                buffer.skip();
                for (unsigned input = 0; input < size; input++) {
                    for (unsigned j = 0; j < size; j++) {
                        expect_true(buffer.isValid(j) == (j != size-input-1));
                    }
                    buffer.put(1);
                }
            }

            test_that("invalidate works") {
                buffer.invalidate(3);
                buffer.invalidate(7);
                for (unsigned j = 0; j < size; j++) {
                    expect_true(!buffer.isValid(j) == (j == 3 || j == 7));
                }
            }
        }
    }
}
