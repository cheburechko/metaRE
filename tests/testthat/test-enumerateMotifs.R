library(testthat)

context("enumerateMotifs")

test_sequences <- c(
    gene1='aaaatgtcaaaa',
    gene2='ccccaaaagggg',
    gene3='ttttggggcccc'
)
test_genes <- names(test_sequences)
k <- 4

test_that("enumerateOligomers with 'genes' output", {
    result <- enumerateOligomers(test_sequences, k, rc=FALSE)
    expect_equal(attr(result, 'geneNames'), names(test_sequences))
    expect_equal(result$`AAAA`, c(1,2))
    expect_equal(result$`CCCC`, c(2,3))
    expect_equal(result$`GGGG`, c(2,3))
    expect_equal(result$`TTTT`, c(3))
    expect_equal(result$`TGTC`, c(1))

    result <- enumerateOligomers(test_sequences, k, rc=TRUE)
    expect_equal(attr(result, 'geneNames'), names(test_sequences))
    expect_equal(result$`AAAA | TTTT`, c(1,2, 3))
    expect_equal(result$`CCCC | GGGG`, c(2,3))
    expect_equal(result$`GACA | TGTC`, c(1))
})

test_that("enumerateOligomers with 'counts' output", {
    result <- enumerateOligomers(test_sequences, k, rc=FALSE, output='counts')
    answers        <- c(3,      2,      2,      1,       1)
    names(answers) <- c('AAAA', 'CCCC', 'GGGG', 'TTTT', 'TGTC')
    expect_equal(result[names(answers)], answers)

    answers        <- c(4,              4,             1)
    names(answers) <- c('AAAA | TTTT', 'CCCC | GGGG', 'GACA | TGTC')
    result <- enumerateOligomers(test_sequences, k, rc=TRUE, output='counts')
    expect_equal(result[names(answers)], answers)
})

test_that("enumerateOligomers with 'positions' output", {
    result <- enumerateOligomers(test_sequences, k, rc=FALSE, output='positions')
    expect_equal(sort(names(result$`AAAA`)), test_genes[1:2])
    expect_equal(result$AAAA$gene1, c(4, 12))
    expect_equal(result$AAAA$gene2, 8)

    expect_equal(sort(names(result$CCCC)), test_genes[2:3])
    expect_equal(result$CCCC$gene2, 4)
    expect_equal(result$CCCC$gene3, 12)

    expect_equal(sort(names(result$GGGG)), test_genes[2:3])
    expect_equal(result$GGGG$gene2, 12)
    expect_equal(result$GGGG$gene3, 8)

    expect_equal(names(result$TTTT), test_genes[3])
    expect_equal(result$TTTT$gene3, 4)

    expect_equal(names(result$TGTC), test_genes[1])
    expect_equal(result$TGTC$gene1, 8)

    result <- enumerateOligomers(test_sequences, k, rc=TRUE, output='positions')
    expect_equal(sort(names(result$`AAAA | TTTT`)), test_genes[1:3])
    expect_equal(result$`AAAA | TTTT`$gene1, c(4, 12))
    expect_equal(result$`AAAA | TTTT`$gene2, 8)
    expect_equal(result$`AAAA | TTTT`$gene3, 4)

    expect_equal(sort(names(result$`CCCC | GGGG`)), test_genes[2:3])
    expect_equal(result$`CCCC | GGGG`$gene2, c(4, 12))
    expect_equal(result$`CCCC | GGGG`$gene3, c(8, 12))

    expect_equal(names(result$`GACA | TGTC`), test_genes[1])
    expect_equal(result$`GACA | TGTC`$gene1, 8)
})

test_that("enumerateOligomers with 'composition' output", {
    result <- enumerateOligomers(test_sequences, k, rc=FALSE, output='composition')
    expect_equal(sort(names(result)), test_genes)

    expect_equal(
        result$gene1,
        c('AAAA', 'AAAT', 'AATG', 'ATGT', 'TGTC', 'GTCA', 'TCAA', 'CAAA', 'AAAA')
    )

    expect_equal(
        result$gene2,
        c('CCCC', 'CCCA', 'CCAA', 'CAAA', 'AAAA', 'AAAG', 'AAGG', 'AGGG', 'GGGG')
    )

    expect_equal(
        result$gene3,
        c('TTTT', 'TTTG', 'TTGG', 'TGGG', 'GGGG', 'GGGC', 'GGCC', 'GCCC', 'CCCC')
    )
})

test_that("enumeratePatterns with 'positions' output", {
    result <- enumeratePatterns(test_sequences, c('AAAA', 'CCCC'), rc=FALSE,
                                output='positions')

    expect_equal(sort(names(result$`AAAA`)), test_genes[1:2])
    expect_equal(result$AAAA$gene1, c(4, 12))
    expect_equal(result$AAAA$gene2, 8)

    expect_equal(sort(names(result$CCCC)), test_genes[2:3])
    expect_equal(result$CCCC$gene2, 4)
    expect_equal(result$CCCC$gene3, 12)

    result <- enumeratePatterns(test_sequences, c('AAAA', 'CCCC'), rc=TRUE,
                                output='positions')
    expect_equal(sort(names(result$AAAA)), test_genes[1:3])
    expect_equal(result$AAAA$gene1, c(4, 12))
    expect_equal(result$AAAA$gene2, 8)
    expect_equal(result$AAAA$gene3, 4)

    expect_equal(sort(names(result$CCCC)), test_genes[2:3])
    expect_equal(result$CCCC$gene2, c(4, 12))
    expect_equal(result$CCCC$gene3, c(8, 12))
})

test_that("enumerateDyadsWithCore with 'positions' output", {
    result <- enumerateDyadsWithCore(test_sequences, k, 'AAAA', -1, 4,
                                           output='positions')
    expect_equal(length(result), 19)
    expect_equal(result$`AAAA_-1_ATGT`$gene1, 7)
    expect_equal(result$`AAAA_0_TGTC`$gene1, 8)
    expect_equal(result$`AAAA_0_GGGG`$gene2, 12)
    expect_equal(result$`CCCC_0_AAAA`$gene3, 8)

    result <- enumerateDyadsWithCore(test_sequences, k, 'TGTC', -1, 4,
                                           fuzzyOrder = TRUE,
                                           output='positions')
    expect_equal(result$`TGTC_0_AAAA`$gene1, c(8,12))

    result <- enumerateDyadsWithCore(test_sequences, k, 'GGGG', -1, 4,
                                           fuzzyOrientation = TRUE,
                                           output='positions')
    expect_equal(result$`AAAA_0_GGGG`$gene2, c(8, 12))
    expect_equal(result$`AAAA_0_GGGG`$gene3, 8)

    result <- enumerateDyadsWithCore(test_sequences, k, 'CCCC', -1, 4,
                                           fuzzySpacer = TRUE,
                                           fuzzyOrientation = TRUE,
                                           output='positions')
    expect_equal(result$`CCCC_-1..4_AAAA`$gene2, c(8,12))
    expect_equal(result$`CCCC_-1..4_AAAA`$gene3, 8)
})
