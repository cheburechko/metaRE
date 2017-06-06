#' @name prepareGEO
#' @title Prepare GEO Samples for differential expression analysis.
#' @description Download given control and treatment samples from NCBI GEO
#' database, convert to log2 scale if necessary.
#' @param control GSM IDs of control samples.
#' @param treatment GSM IDs of treatment samples.
#' @param isLog2 logical flag which indicates whether the original data is in
#' log2 scale.
#' @param id_col column name for row ids
#' @param data_col column name for expression values
#' @return a data frame with a column for each sample with log2 expression
#' levels for each gene.
#' @examples
#' prepareGEO(
#'     c("GSM469656", "GSM469657", "GSM469658"),
#'     c("GSM469650", "GSM469651", "GSM469652"),
#'     F
#' )
#' @export
#' @importFrom Biobase assayData
#' @importFrom GEOquery getGEO GSMList Table
#' @importFrom stats setNames
prepareGEO <- function(control, treatment, isLog2, id_col='ID_REF', data_col='VALUE') {
    samples <- c(control, treatment)

    data <- NULL
    for (sample in samples) {
        if (is.null(data)) {
            data <- Table(getGEO(sample))
            rownames(data) <- as.character(data[, id_col])
            data[, sample] <- data[, data_col]
            data <- data[, sample, drop=F]
        } else {
            buf <- Table(getGEO(sample))
            data[as.character(buf[, id_col]), sample] <- buf[, data_col]
        }
    }

    if (!isLog2) {
        data <- log2(data)
    }

    return(data)
}
