library(testthat)

context('MassContingencyTable')

test_that("MassContingencyTable", {
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
    test <- calculateMassContingencyTablePvalues(gcs, gcm)
    expect_equal(result, test)

    expect_error(calculateMassContingencyTablePvalues(geneStructure, gcm))
    expect_error(calculateMassContingencyTablePvalues(gcs, mat))
})
