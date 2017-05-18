#' @name GeneClassificationSparse
#' @title Gene Classification With Sparse Structure
#' @description Wrap a list of integer vectors, which indicate the genes in
#' which the motifs are present for use in
#' \code{\link{calculateMassContingencyTablePvalues}}
#' @param x a named list of integer vectors, each vector is named after a motif,
#' integers describe indice for genes in which the motif was detected.
#' @param geneNames character vector with gene names which describes indices in
#' \code{x}.
#' @param gcs an object of 'GeneClassificationSparse' class
#' @details
#' \code{GeneClassificationSparse} adds \code{geneNames} as attribute to
#' \code{x} and adds \code{'GeneClassificationSparse'} to classes of \code{x}
#' @return
#' \code{GeneClassificationMatrix} returns new \code{GeneClassificationSparse}
#' object
#'
#' \code{geneNames} returns geneNames from \code{GeneClassificationSparse}
#' object
#' @examples
#'  x <- list(
#'     elem1=c(1,2,3),
#'     elem2=c(2, 5),
#'     elem3=c(10, 6)
#' )
#' genes <- paste0('gene', 1:10)
#' gcs <- GeneClassificationSparse(x, genes)
#'
#' gcs
#' geneNames(gcs)
#' @export
GeneClassificationSparse <- function(x, geneNames) {
    if (!inherits(x, 'list') || !all(sapply(x, is.numeric))) {
        stop("x must be the list of integer vectors")
    }
    if (!is.character(geneNames)) {
        stop("geneNames must be a character vector")
    }
    if (min(sapply(x, min)) < 1 || max(sapply(x, max)) > length(geneNames)) {
        stop("integers in x must be in range [1, length(geneNames)]")
    }
    attr(x, 'geneNames') <- geneNames
    class(x) <- c('GeneClassificationSparse', class(x))
    x
}

#' @rdname GeneClassificationSparse
#' @export
geneNames <- function(gcs) {
    if (!inherits(gcs, 'GeneClassificationSparse')) {
        stop("x must have class ''GeneClassificationSparse'")
    }
    attr(gcs, 'geneNames')
}

#' @name GeneClassificationMatrix
#' @title Gene Classification With Dense Structure
#' @description Wrap a logical matrix which describes which genes were
#' differentially expressed in the give experiments for use in
#' \code{\link{calculateMassContingencyTablePvalues}}
#' @param mat a logical matrix with rows corresponding to genes and columns
#' corresponding to experiments.
#' @param x a \code{GeneClassificationMatrix} object
#' @details \code{GeneClassificationMatrix} adds \code{geneCounts} attribute to
#' \code{mat} which describes column sums of the matrix and adds
#' \code{'GeneClassificationMatrix'} to \code{mat} classes.
#' @return
#' \code{GeneClassificationMatrix} returns new \code{GeneClassificationMatrix}
#' object
#'
#'
#' \code{geneCounts} returns geneCounts from \code{GeneClassificationMatrix}
#' @examples
#' data <- matrix(runif(30) < 0.5, 6, 5)
#' rownames(data) <- paste0('gene', 1:6)
#' colnames(data) <- paste0('exp', 1:5)
#' gcm <- GeneClassificationMatrix(data)
#'
#' gcm
#' geneCounts(gcm)
#' @export
GeneClassificationMatrix <- function(mat) {
    if (!inherits(mat, c('matrix', 'GeneClassificationMatrix')) || !is.logical(mat)) {
        stop('mat must be a logical matrix')
    }
    attr(mat, 'geneCounts') <- colSums(mat)
    class(mat) <- c('GeneClassificationMatrix', class(mat))
    mat
}

#' @rdname GeneClassificationMatrix
#' @export
geneCounts <- function(x) {
    if (!inherits(x, 'GeneClassificationMatrix')) {
        stop("x must have class ''GeneClassificationMatrix'")
    }
    attr(x, 'geneCounts')
}
