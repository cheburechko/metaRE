#' @export
#' @importFrom Biobase assayData
#' @importFrom GEOquery getGEO GSMList
#' @importFrom stats setNames
prepareGEO <- function(control, treatment, isLog2, filename=NULL, GEO=NULL) {
    if (is.null(filename) && is.null(GEO)) {
        stop("Either filename or GEO ID must be provided")
    }
    #eset <- getGEO(GEO=GEO, filename=filename)
    samples <- c(control, treatment)

    data <- NULL
    for (sample in samples) {
        if (is.null(data)) {
            data <- Table(getGEO(sample))
            rownames(data) <- data$ID_REF
            data[, sample] <- data$VALUE
            data <- data[, sample, drop=F]
        } else {
            data[, sample] <- Table(getGEO(sample))$VALUE
        }
    }

    if (!isLog2) {
        data <- log2(data)
    }

    return(data)
}
