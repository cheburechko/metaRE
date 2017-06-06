#' @title Permutation test on Meta P-Values of Association
#' @name permutationTest
#' @description Calculate permutation test in order to evaluate robustness of
#' the meta p-values obtained in \code{\link{testRegulationHypotheses}}.
#' @param hypothesesClasses An object of \code{\link{GeneClassificationSparse}}
#' class, usually describes sets of genes in which motifs are present
#' @param annotationClasses An object of \code{\link{GeneClassificationMatrix}}
#' class, usually describes sets of genes which were differentially expressed in
#' experiments
#' @param n number of permutations
#' @param alternative indicates the alternative hypothesis and must be one of
#' "two.sided", "greater" or "less"
#' @param outfile path to file for writing the preliminary results of permuation
#' test
#' @param pvaluePreFilter threshold on permutation p-value: if the
#' hypotheses gets p-value above the threshold, it is excluded from the
#' analysis.
#' @param perRun number of permuations for a single run (see Details)
#'
#' @details Each gene has a specific hypotheses profile (described by
#' \code{hypothesesClasses}) and a specific annotation profile (described by
#' \code{annotationClasses}). One can see a gene as a connection between two
#' profiles. The permutation test breaks this connection and assigns a random
#' hypotheses profile to each annotation profile without replacement. (E.g. if
#' hypotheses are regulatory elements in gene promoters and annotations are
#' significant differential expressions of genes in a set of experiments,
#' one can say that permutation test assigns a random promoter to a gene). This
#' way test permutation test applies random permutations while keeping the
#' profiles realistic from biological standpoint.
#'
#' The permutation p-value is calculated as \code{perm.p.value = (LEQ+1)/(n+1)},
#' where \code{LEQ} is the number of permutations where meta p-values were less
#' or equal to the real meta p-values, \code{n} is the number of permutations.
#'
#' The test automatically writes the intermediate results to \code{outfile} and
#' the previous intermediate results to \code{outfile} with '.last' suffix
#' and excludes hypotheses that will have the permutation p-value above
#' \code{pvaluePreFilter} from analysis after every \code{perRun} permutations.
#'
#' \code{permuationTest} uses \code{\link{foreach}} and will run in parallel
#' if a parallel foreach backend is registered (e.g.
#' \code{\link[doParallel]{registerDoParallel}}). The calculations will be
#' performed in batches of size \code{perRun}, which will be split equally
#' between workers.
#'
#' @return A \code{\link{data.frame}} with the following columns:
#' \describe{
#' \item{Hypothesis}{hypothesis name taken from \code{hypothesesClasses}}
#' \item{LEQ}{number of permutations where meta p-value was less or equal to the
#' real meta p-value.}
#' \item{Permuations}{total number of permutations}
#' \item{Meta.P.Value}{raw meta p-values}
#' \item{Permutation.P.Value}{raw permuation p-values}
#' }
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
#' output <- tempfile()
#' test <- permutationTest(gcs, gcm, 1000, outfile=output, perRun=100)
#' print(test)
#'
#' \dontrun{
#' ## Using doParallel
#'
#' library(doParallel)
#' registerDoParallel()
#' test <- permutationTest(gcs, gcm, 1000, outfile=output, perRun=500)
#' print(test)
#' }
#'
#' @importFrom futile.logger flog.info
#' @importFrom foreach getDoParWorkers %dopar% foreach
#' @importFrom utils write.csv
#' @export
permutationTest <- function(
    hypothesesClasses, annotationClasses, n,
    alternative=c('greater', 'less', 'two.sided'),
    outfile='./perm_test.csv',
    pvaluePreFilter=NULL,
    perRun=n
) {
    alternative <- match.arg(alternative)
    preparedData <- .preparePermutaionData(hypothesesClasses, annotationClasses,
                                           alternative, outfile)
    done <- 0
    run <- 0
    while (done < n) {
        run <- run + 1
        flog.info("Starting run %d", run)
        thisRun <- min(perRun, n-done)
        preparedData$result <- preparedData$result + .permutationTestRun(
            preparedData$hypothesesClasses, preparedData$annotationClasses,
            preparedData$realMetaPValues, alternative, thisRun
        )
        done <- done + thisRun
        preparedData <- .dropExtraElements(preparedData, n, pvaluePreFilter)

        if (file.exists(preparedData$cur_file)) {
            file.rename(preparedData$cur_file, preparedData$last_file)
        }

        permPValue <- .permutationPvalue(preparedData$result, done)
        df <- data.frame(
            Hypothesis=names(preparedData$hypothesesClasses),
            LEQ=preparedData$result,
            Permutations=done,
            Meta.P.Value=preparedData$realMetaPValues,
            Permutation.P.Value=permPValue,
            row.names = names(preparedData$hypothesesClasses),
            stringsAsFactors=FALSE
        )
        df <-df[order(df$Permutation.P.Value, df$Meta.P.Value), ]
        write.csv(df, preparedData$cur_file, row.names = FALSE)
        flog.info("Hypotheses left: %d", nrow(df))
    }
    return(df)
}

.permutationPvalue <- function(leq, total) (leq+1)/(total+1)

.preparePermutaionData <- function(hypothesesClasses, annotationClasses,
                                   alternative, outfile) {
    metaAnalysis <- calcMetaAssociation(
        calculateMassContingencyTablePvalues(hypothesesClasses, annotationClasses, alternative)
    )[names(hypothesesClasses), ]
    annotationClasses<-.reorderAnnotationGenes(hypothesesClasses, annotationClasses)
    list(
        hypothesesClasses=hypothesesClasses,
        annotationClasses=annotationClasses,
        realMetaPValues=metaAnalysis$Meta.P.Value,
        last_file=paste0(outfile, '.last'),
        cur_file=outfile,
        result=rep(0, length(hypothesesClasses))
    )
}

.dropExtraElements <- function(preparedData, n, threshold=NULL) {
    output <- preparedData
    if (is.null(threshold)) {
        return(output)
    }
    index <- .permutationPvalue(preparedData$result, n) <= threshold
    output$hypothesesClasses <- preparedData$hypothesesClasses[index]
    output$result <- preparedData$result[index]
    output$realMetaPValues <- preparedData$realMetaPValues[index]

    return(output)
}

.permutationTestRun <- function(hypothesesClasses, annotationClasses,
                                realMetaPValues, alternative, n)
{
    cores <- getDoParWorkers()
    perCore <- ceiling(n/cores)
    tasks <- lapply(1:cores, function(core) ((core-1)*perCore+1):min(n, core*perCore))
    foreach(
        task=tasks,
        .export=c('hypothesesClasses', 'annotationClasses', 'realMetaPValues',
                  'alternative'),
        .inorder = FALSE,
        .combine=`+`
    ) %dopar% {
        m <- rep(0, length(hypothesesClasses))
        for (i in task) {
            pvalues <- massFisherTest(
                annotationClasses[sample(1:nrow(annotationClasses)), ],
                geneCounts(annotationClasses),
                hypothesesClasses, alternative
            )
            simPvalues <- bulkSumlog(pvalues)
            m <- m + (simPvalues <= realMetaPValues)
        }
        return(m)
    }
}
