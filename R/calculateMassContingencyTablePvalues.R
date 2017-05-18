.reorderAnnotationGenes <- function(hypothesesClasses, annotationClasses) {
    GeneClassificationMatrix(annotationClasses[geneNames(hypothesesClasses), ])
}

#' @name MassContingencyTable
#' @title Calculate P-Values For Many Contingency Tables
#' @description Calculate p-values for a set of contingency tables, defined by
#' all possible combinations between two different sets of object
#' classifications using the Fisher's exact test.
#' @param hypothesesClasses An object of \code{\link{GeneClassificationSparse}}
#' class, usually describes sets of genes in which motifs are present
#' @param annotationClasses An object of \code{\link{GeneClassificationMatrix}}
#' class, usually describes sets of genes which were differentially expressed in
#' experiments
#' @param alternative indicates the alternative hypothesis and must be one of
#' "two.sided", "greater" or "less"
#' @details \code{geneNames(hypothesesClasses)} and
#' \code{rownames(annotationClasses)} must describe the same set of genes. Gene
#' order though can be different, the function reorders the genes itself.
#' @return A float matrix with p-values. Columns correspond to columns in
#' \code{annotationClasses}, rows correspond to items in
#' \code{hypothesesClasses}.
#' @examples
#' elements <- 5
#' genes <- 200
#' experiments <- 5
#'
#' geneNames <- paste0('gene', 1:genes)
#' expNames <- paste0('exp', 1:experiments)
#' elemNames <- paste0('elem', as.character(1:elements))
#'
#' gcm <- GeneClassificationMatrix(
#'     matrix(runif(genes*experiments)<0.05, genes, experiments,
#'     dimnames=list(geneNames, expNames))
#' )
#'
#' gcs <- GeneClassificationSparse(
#'     setNames(lapply(1:elements, function(x) {
#'         sample(1:genes, as.integer(runif(1, max=genes)))
#'     }), elemNames),
#'     geneNames
#' )
#'
#' test <- calculateMassContingencyTablePvalues(gcs, gcm)
#' print(test)
#' @export
calculateMassContingencyTablePvalues <- function(
    hypothesesClasses, annotationClasses,
    alternative=c('greater', 'less', 'two.sided')
) {
    if (!inherits(hypothesesClasses, 'GeneClassificationSparse')) {
        stop("hypothesesClasses must have 'GeneClassificationSparse' class")
    }
    if (!inherits(annotationClasses, 'GeneClassificationMatrix')) {
        stop("annotationClasses must have 'GeneClassifcationMatrix' class")
    }

    annotationClasses <- .reorderAnnotationGenes(hypothesesClasses, annotationClasses)
    alternative <- match.arg(alternative)

    result <- massFisherTest(annotationClasses, geneCounts(annotationClasses),
                             hypothesesClasses, alternative)

    dimnames(result) <- list(names(hypothesesClasses), colnames(annotationClasses))

    result
}
