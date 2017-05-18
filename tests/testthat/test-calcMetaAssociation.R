library(testthat)

context("calcMetaAssociation")

test_that("bulkSumlog", {
    cols <- 10
    rows <- 10
    data <- matrix(runif(cols*rows), rows, cols)
    expect_equal(
        bulkSumlog(data),
        apply(data, 1, function(x) {
            chisq <- -2*sum(log(x))
            pchisq(chisq, 2*length(x), lower.tail = F)
        })
    )
})

test_that("bulkSumlog with bad values", {
    cols <- 10
    rows <- 10
    data <- matrix(runif(cols*rows), rows, cols)
    data[sample(1:(cols*rows), 5)] <- 0
    data[sample(1:(cols*rows), 5)] <- 2
    expect_warning(result <- bulkSumlog(data), 'P-values must be')

    data[data == 0] <- .Machine$double.xmin / 2**(.Machine$double.digits-1)
    data[data > 1] <- 1

    expect_equal(
        result,
        apply(data, 1, function(x) {
            chisq <- -2*sum(log(x))
            pchisq(chisq, 2*length(x), lower.tail = F)
        })
    )
})

test_that("bulkSumlog with single column", {
    rows <- 10
    data <- runif(rows)
    expect_equal(
        bulkSumlog(data),
        data
    )
})

test_that("calcMetaAssociaiton", {
    elements <- 1000
    experiments <- 5
    elemNames <- paste0('elem', 1:elements)
    expNames <- paste0('exp', 1:experiments)
    assocTable <- matrix(runif(elements*experiments), elements, experiments,
                         dimnames=list(elemNames, expNames))

    result <- calcMetaAssociation(assocTable, 'bonferroni')
    expect_is(result, 'MetaAssociationTable')
    expect_equal(attr(result, 'adjustMethod'), 'bonferroni')
    expect_equal(attr(result, 'experimentCount'), experiments)
    expect_equal(length(intersect(rownames(result), elemNames)), elements)
    p_values <- sort(bulkSumlog(assocTable))
    names(p_values) <- NULL
    expect_equal(p_values, result$Meta.P.Value)
    expect_equal(p.adjust(p_values, 'bonferroni'), result$Adj.Meta.P.Value)
})

test_that("testRegulationHypotheses", {
    elements <- 3
    genes <- 200
    experiments <- 3
    genesWithElem <- 70
    genesWithoutElem <- genes - genesWithElem
    degs <- c(50, 100, 150)
    degsWithElem <- c(20, 30, 40)

    geneNames <- paste0('gene', 1:genes)
    expNames <- paste0('exp', 1:experiments)
    names(degs) <- expNames
    mat <- matrix(FALSE, genes, experiments, dimnames=list(geneNames, expNames))
    for (exp in 1:experiments) {
        mat[, exp] <- c(rep(TRUE, degs[exp]), rep(FALSE, genes-degs[exp]))
    }
    gcm <- GeneClassificationMatrix(mat)

    elemNames <- paste0('elem', as.character(1:elements))
    geneStructure <- list()
    for (elem in 1:elements) {
        geneStructure[[elemNames[elem]]] <- c(
            1:degsWithElem[elem], (degsWithElem[elem]+genesWithoutElem+1):genes
        )
    }
    gcs <- GeneClassificationSparse(geneStructure, geneNames)
    result <- matrix(1, elements, experiments,
                     dimnames = list(elemNames, expNames))

    for (exp in 1:experiments) {
        for (elem in 1:elements) {
            result[elem, exp] <- fisher.test(matrix(c(
                degsWithElem[elem], genesWithElem-degsWithElem[elem],
                degs[exp]-degsWithElem[elem],
                genes-genesWithElem-degs[exp]+degsWithElem[elem]
            ), 2, 2), alternative = 'greater')$p.value
        }
    }
    metaResult <- bulkSumlog(result)

    test <- testRegulationHypotheses(gcs, gcm, adjust='bonferroni')
    metaResult <- setNames(metaResult[rownames(test)], NULL)

    expect_is(test, 'MetaAssociationTable')
    expect_equal(attr(test, 'adjustMethod'), 'bonferroni')
    expect_equal(attr(test, 'experimentCount'), experiments)

    expect_equal(length(intersect(rownames(test), elemNames)), elements)

    expect_equal(metaResult, test$Meta.P.Value)
    expect_equal(p.adjust(metaResult, 'bonferroni'), test$Adj.Meta.P.Value)
})
