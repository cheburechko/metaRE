#' @name exampleData
#' @title Example Datasets
#' @description This is the example data that is used in vignette:
#' \code{vignette("quickstart", package = "metaRE")}
#'
#' \describe{
#'     \item{MA_AT_auxin}{List of NCBI GEO microarray datasets for auxin
#'     treatment of \emph{A. thaliana}}
#'     \item{AT_AGI}{List of AGI identifiers of genes used in analysis}
#' }
"MA_AT_auxin"

MA_AT_auxin <- list(
    list(
        control=c("GSM1030559", "GSM1030575", "GSM1030591"),
        treatment=c("GSM1030567", "GSM1030583", "GSM1030599"),
        log2=T, name="GSE42007_0.5h"
    ),
    list(
        control=c("GSM1030560", "GSM1030576", "GSM1030592"),
        treatment=c("GSM1030568", "GSM1030584", "GSM1030600"),
        log2=T, name="GSE42007_1h"
    ),
    list(
        control=c("GSM9571", "GSM9572", "GSM9573"),
        treatment=c("GSM9574", "GSM9575", "GSM9576"),
        log2=T, name="GSE627"
    ),
    list(
        control=c("GSM25858", "GSM25859", "GSM25860", "GSM25861"),
        treatment=c("GSM25862", "GSM25863", "GSM25864", "GSM25865"),
        log2=T, name="GDS1044"
    ),
    list(
        control=c("GSM18228", "GSM18230", "GSM18232", "GSM18229", "GSM18231", "GSM18233"),
        treatment=c("GSM18290", "GSM18291", "GSM18292", "GSM18293"),
        log2=F, name="GDS672_0.1uM_1h"
    ),
    list(
        control=c("GSM18228", "GSM18230", "GSM18232", "GSM18229", "GSM18231", "GSM18233"),
        treatment=c("GSM18298", "GSM18299", "GSM18300", "GSM18301"),
        log2=F, name="GDS672_1uM_1h"
    ),
    list(
        control=c("GSM13430", "GSM13432"),
        treatment=c("GSM13433", "GSM13434"),
        log2=F, name="GDS744"
    ),
    list(
        control=c("GSM1053036", "GSM1053037", "GSM1053038"),
        treatment=c("GSM1053030", "GSM1053031", "GSM1053032"),
        log2=T, name="GSE42896"
    ),
    list(
        control=c("GSM469656", "GSM469657", "GSM469658"),
        treatment=c("GSM469653", "GSM469654", "GSM469655"),
        log2=F, name="GSE18975_0.5h"
    ),
    list(
        control=c("GSM469656", "GSM469657", "GSM469658"),
        treatment=c("GSM469650", "GSM469651", "GSM469652"),
        log2=F, name="GSE18975_1h"
    )
)

#' @rdname exampleData
"AT_AGI"
