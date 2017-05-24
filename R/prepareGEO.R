#' @export
#' @importFrom Biobase assayData
#' @importFrom GEOquery getGEO GSMList
#' @importFrom stats setNames
prepareGEO <- function(control, treatment, isLog2, filename=NULL, GEO=NULL) {
    if (is.null(filename) && is.null(GEO)) {
        stop("Either filename or GEO ID must be provided")
    }

    eset <- getGEO(GEO=GEO, filename=filename)
    samples <- c(control, treatment)

    if (class(eset) == "GDS") {
        data <- eset@dataTable@table[, samples]
        rownames(data) <- eset@dataTable@table$ID_REF
    } else if (class(eset) == "GSE") {
        gsmlist <- setNames(lapply(
            GSMList(eset)[samples],
            function(x) x@dataTable@table
        ), samples)
        data <- data.frame(row.names=as.character(gsmlist[[1]]$ID_REF))
        for (name in names(gsmlist)) {
            data[as.character(gsmlist[[name]]$ID_REF), name] <-
                gsmlist[[name]]$VALUE
        }
    } else {
        if (is.list(eset)) {
            eset <- eset[[1]]
        }
        data <- as.data.frame(assayData(eset)$exprs[, samples])
    }

    if (!isLog2) {
        data <- log2(data)
    }

    return(data)
}
