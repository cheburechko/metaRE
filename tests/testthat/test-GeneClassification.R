library(testthat)

context('GeneClassifcation')

test_that("GeneClassifcationSparse", {
    x <- list(
        elem1=c(1,2,3),
        elem2=c(2, 5),
        elem3=c(10, 6)
    )
    genes <- paste0('gene', 1:10)
    gcs <- GeneClassificationSparse(x, genes)

    expect_equal(genes, geneNames(gcs))
    expect_equivalent(x, gcs)

    expect_error(
        GeneClassificationSparse(c(10, 20, 30), paste0('gene', 1:30)),
        'list of integer vectors'
    )
    expect_error(
        GeneClassificationSparse(list(elem1='1'), genes),
        'list of integer vectors'
    )
    expect_error(
        GeneClassificationSparse(list(elem1=c(10,20,40)), genes),
        'must be in range'
    )
    expect_error(
        GeneClassificationSparse(x, 1:10),
        "geneNames must be a character vector"
    )
    expect_error(geneNames(x))
})

test_that("GeneClassificationMatrix", {
    data <- matrix(runif(30) < 0.5, 6, 5)
    rownames(data) <- paste0('gene', 1:6)
    colnames(data) <- paste0('exp', 1:5)
    gcm <- GeneClassificationMatrix(data)

    expect_equal(geneCounts(gcm), colSums(data))
    expect_equivalent(gcm[1:length(data)], data[1:length(data)])

    expect_error(GeneClassificationMatrix(matrix(runif(30), 6, 5)))
    expect_error(geneCounts(data))
})
