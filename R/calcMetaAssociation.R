#' @title Combine p-values by the sum of logs method (bulk)
#' @description Combine p-values by the sum of logs method, also known as
#' Fisher's method
#' @param m A matrix of p-values (0 < p <= 1)
#' @details Implementation of \code{\link[metap]{sumlog}} for matrices. Combines
#' p-values in rows of the matrix.
#' @return A float vector of meta p-values corresponding to each row of the
#' initial matrix.
#' @examples
#' data <- matrix(runif(100), 10, 10)
#' bulkSumlog(data)
#' @export
bulkSumlog = function(m) {
    if (is.vector(m)) {
        return(m)
    }
    if (any(m == 0)) {
        warning("P-values must be > 0")
        m[m == 0] <- .Machine$double.xmin / 2**(.Machine$double.digits-1)
    }
    if (any(m > 1)) {
        warning("P-values must be <= 1")
        m[m > 1] <- 1
    }
    chisq <- -2 * rowSums(log(m))
    df <- 2 * ncol(m)
    pchisq(chisq, df, lower.tail = FALSE)
}

#' @title Calculate Meta P-Values Of Association
#' @description Calculate meta p-values of association for every motif across
#' all experiments with Fisher's method
#' @param assocTable a double matrix of p-values, result of
#' \code{calculateMassContingencyTablePvalues}
#' @param adjust multiple testing correction method
#' (see \code{\link{p.adjust.methods}})
#' @param threshold cutoff value for adjusted p-value
#' @return A \code{data.frame} with class \code{'MetaAssociationTable'}. It has
#' the following columns:
#' \describe{
#' \item{Meta.P.Value}{raw meta p-values}
#' \item{Adj.Meta.P.Value}{meta p-values adjusted with one of the multiple
#' testing correction methods}
#' }
#'
#' It also has the following attributes
#' \describe{
#' \item{adjustMethod}{method used for multiple testing correciton}
#' \item{experimentCount}{number of experiments used in analysis}
#' }
#' @export
#' @examples
#' elements <- 10
#' experiments <- 10
#' elemNames <- paste0('elem', 1:elements)
#' expNames <- paste0('exp', 1:experiments)
#' data <- matrix(runif(elements*experiments), elements, experiments,
#'                dimnames=list(elemNames, expNames))
#' calcMetaAssociation(data, 'fdr')
#' @importFrom stats p.adjust pchisq
calcMetaAssociation <- function(assocTable, adjust='fdr', threshold=1.) {
    df <- data.frame(Meta.P.Value=rep(1, nrow(assocTable)),
                     Adj.Meta.P.Value=1,
                     row.names=rownames(assocTable),
                     stringsAsFactors=FALSE)

    df$Meta.P.Value <- bulkSumlog(assocTable)
    df$Adj.Meta.P.Value <- p.adjust(df$Meta.P.Value, adjust)
    df <- df[df$Adj.Meta.P.Value <= threshold, ]
    df <- df[order(df$Adj.Meta.P.Value, df$Meta.P.Value), ]

    attr(df, 'adjustMethod') <- adjust
    attr(df, 'experimentCount') <- ncol(assocTable)
    class(df) <- c('MetaAssociationTable', class(df))
    df
}

#' @title Test hypotheses for association with gene regulation
#' @name testRegulationHypotheses
#' @description Test given set of hypotheses for association with gene
#' regulation performing meta-analysis over many experiments.
#' @param hypothesesClasses An object of \code{\link{GeneClassificationSparse}}
#' class, usually describes sets of genes in which motifs are present
#' @param annotationClasses An object of \code{\link{GeneClassificationMatrix}}
#' class, usually describes sets of genes which were differentially expressed in
#' experiments
#' @param alternative indicates the alternative hypothesis and must be one of
#' "two.sided", "greater" or "less"
#' @param adjust multiple testing correction method
#' @param threshold cutoff value for adjusted p-value
#' (see \code{\link{p.adjust.methods}})
#' @details This is basically a shorthand for
#' \code{\link{calculateMassContingencyTablePvalues}} followed by
#' \code{\link{calcMetaAssociation}}
#' @return Same as \code{\link{calcMetaAssociation}}.
#' @seealso \code{\link{calculateMassContingencyTablePvalues}},
#' \code{\link{calcMetaAssociation}}
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
#' test <- testRegulationHypotheses(gcs, gcm)
#' print(test)
#' @export
testRegulationHypotheses <- function(
    hypothesesClasses, annotationClasses,
    alternative=c('greater', 'less', 'two.sided'),
    adjust='fdr', threshold=1.
) {
    assocTable <- calculateMassContingencyTablePvalues(
        hypothesesClasses, annotationClasses, alternative
    )
    calcMetaAssociation(assocTable, adjust, threshold)
}
