library(testthat)

context("permutationTest")

elements <- paste0('elem', 1:10)
genes <- paste0('gene', 1:100)
experiments <- paste0('exp', 1:10)
percentage <- 0.05

gcm <- matrix(
    runif(length(genes)*length(experiments)) < percentage,
    length(genes), length(experiments), dimnames=list(genes, experiments)
)
gcm[1:10, ] <- TRUE
gcm <- GeneClassificationMatrix(gcm)

gcs <- matrix(
    runif(length(genes)*length(elements)) < percentage,
    c(length(genes), length(elements))
)
gcs[1:10,  ] <- FALSE
gcs[1:10, 1] <- TRUE

element_indices <- setNames(apply(gcs, 2, which), elements)
gcs <- GeneClassificationSparse(element_indices, genes)

metaPValue <- setNames(bulkSumlog(
    calculateMassContingencyTablePvalues(gcs, gcm)
), NULL)

permutations <- 1000

library(futile.logger)
flog.threshold(WARN)

test_that("permutation test works", {
    tempf <- tempfile(fileext = ".csv")
    foreach::registerDoSEQ()

    result <- permutationTest(gcs, gcm, n=permutations, outfile=tempf)[elements, ]

    expect_equal(nrow(result), length(elements))
    expect_equal(result$Meta.P.Value, metaPValue)
    expect_equal((result$LEQ+1)/(result$Permutations+1), result$Permutation.P.Value)

    file.remove(tempf)
})

test_that("permutation test works with prefilter", {
    tempf <- tempfile(fileext = ".csv")
    foreach::registerDoSEQ()

    result <- permutationTest(
        gcs, gcm, n=permutations, outfile=tempf, pvaluePreFilter=0.05
    )

    expect_equal(rownames(result), elements[1])
    expect_equal(result$Meta.P.Value, metaPValue[1])
    expect_equal((result$LEQ+1)/(permutations+1), result$Permutation.P.Value)

    file.remove(tempf)
})

test_that("permutation test works several runs", {
    tempf <- tempfile(fileext = ".csv")
    foreach::registerDoSEQ()

    result <- permutationTest(
        gcs, gcm, n=permutations, outfile=tempf, perRun=permutations/2
    )

    last_file <- paste0(tempf, '.last')
    last_data <- read.csv(last_file, stringsAsFactors = F)
    expect_equal(last_data$Permutations, rep(permutations/2, length(elements)))
    expect_equal(sort(last_data$Hypothesis), sort(elements))

    latest_data <- read.csv(tempf, stringsAsFactors = F)
    rownames(latest_data) <- latest_data$Hypothesis
    expect_equal(latest_data, result)

    result <- result[elements, ]

    expect_equal(nrow(result), length(elements))
    expect_equal(result$Meta.P.Value, metaPValue)
    expect_equal((result$LEQ+1)/(permutations+1), result$Permutation.P.Value)
    file.remove(tempf)
    file.remove(last_file)
})

test_that("parallel permutation test works", {
    skip_if_not_installed('doParallel')
    doParallel::registerDoParallel(4)
    tempf <- tempfile(fileext = ".csv")

    result <- permutationTest(gcs, gcm, n=permutations, outfile=tempf)[elements, ]
    expect_equal(result$Meta.P.Value, metaPValue)
    expect_equal(nrow(result), length(elements))
    expect_equal((result$LEQ+1)/(permutations+1), result$Permutation.P.Value)

    foreach::registerDoSEQ()
    file.remove(tempf)
})

flog.threshold(INFO)
