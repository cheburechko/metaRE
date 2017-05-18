library(testthat)

context("contTableGenerator")

test_that("massFisherTest works like fisher.test", {
    samples <- 1000
    maxValues <- 1000
    n1 <- ceiling(runif(samples, min=1, max=maxValues))
    n2 <- ceiling(runif(samples, min=1, max=maxValues))
    eff1 <- sapply(n1, function(x) ceiling(runif(1, 0, x)))
    eff2 <- sapply(n2, function(x) ceiling(runif(1, 0, x)))

    for (alternative in c('two.sided', 'less', 'greater')) {
        trueValues <- sapply(1:samples, function(i) {
            fisher.test(matrix(c(eff1[i], eff2[i], n1[i]-eff1[i], n2[i]-eff2[i]), 2, 2), alternative = alternative)$p.value
        })
        testValues <- quickFisherTest(eff1, n1, eff2, n2, alternative)
        expect_equal(trueValues, testValues)

        expect_equal(testValues, trueValues)
    }

})

test_that("massFisherTest works like fisher.test in specific case 1", {
    trueValue <- fisher.test(matrix(c(141, 71, 229-141, 114-71), 2, 2), alternative = "two.sided")$p.value
    testValue <- quickFisherTest(141, 229, 71, 114, "two.sided")
    expect_equal(trueValue, testValue)
})

