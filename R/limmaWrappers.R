#' @importFrom limma topTable
.processFit = function(fit, classes, adjust='none') {
   df <- topTable(fit, number=Inf, adjust.method=adjust)

   result <- data.frame(row.names=rownames(df))
   for (name in names(classes)) {
       result[ ,name] <- classes[[name]](df)
   }
   return(result)
}

#' @name preprocessDEGs
#' @title Extract DEGs from microarrays and RNA-Seq
#' @param df data frame where column names are samples, row names are gene
#' names, data is either expression levels or rna counts of genes in given samples
#' @param treatment character vector, names of treamtment samples
#' @param control character vector, names of control samples
#' @param classes named list of functions, each function must take one argument
#' - a data frame which is returned by \code{\link[limma]{topTable}} function from
#' \code{\link{limma}} package, and must return a logical vector - genes that
#' are considered as DEGs.
#' @param adjust multiple correction method (see \code{\link{p.adjust}})
#' @param dataList named list of experiments where each item is a list with the
#' following items:
#' \describe{
#'     \item{data}{same as \code{df} parameter}
#'     \item{treatment}{same as \code{treament} parameter}
#'     \item{control}{same as \code{control} parameter}
#'     \item{type}{either 'MA' or 'RNA', describes the type of the dataset}
#' }
#' @description Given microarray or RNA-Seq data extract DEGs given user-defined
#' thresholds on p-values and fold changes.
#' @details The search is perfromed using \code{\link{limma}} and
#' \code{\link{edgeR}} packages. User should define the thresholds on DEGs using
#' \code{classes} parameters. The following columns as returned by
#' \code{\link[limma]{topTable}} are suggested to be used for filtering:
#' \describe{
#'     \item{logFC}{log2-fold-change}
#'     \item{P.Value}{raw p-value}
#'     \item{adj.P.Value}{adjusted p-value with \code{adjust} method}
#' }
#' @return \code{processMicroarray} and \code{processRNACounts} return data
#' frame with a column for each class of DEGs as defined in \code{classes}
#' parameter. Each column is a logical vector which tells if a gene (row)
#' belongs to this DEG class.
#'
#' \code{preprocessGeneExpressionData} returns a named list of
#' \code{\link{GeneClassificationMatrix}}, each matrix corresponds to an item in
#' \code{classes}
#' @examples
#' # The following is taken from limma::lmFit examples
#' # Simulate gene expression data for 100 probes and 6 microarrays
#' # Microarray are in two groups
#' # First two probes are differentially expressed in second group
#' # Std deviations vary between genes with prior df=4
#' sd <- 0.3*sqrt(4/rchisq(100,df=4))
#' y <- matrix(rnorm(100*6,sd=sd),100,6)
#' rownames(y) <- paste("Gene",1:100)
#' colnames(y) <- paste("Experiment", 1:6)
#' control <- colnames(y)[1:3]
#' treatment <- colnames(y)[4:6]
#' y[1:2,4:6] <- y[1:2,4:6] + 2
#' ma_data <- y
#'
#' rna_data <- matrix(as.integer(exp(y)*100), 100, 6)
#' rownames(rna_data) <- paste("Gene",1:100)
#' colnames(rna_data) <- paste("Experiment", 1:6)
#'
#' dataList <- list(
#'    exp1=list(data=ma_data, control=control, treatment=treatment, type='MA'),
#'    exp2=list(data=rna_data, control=control, treatment=treatment, type='RNA')
#' )
#'
#' classes <- list(
#'     up=function(df) df$logFC > 0 & df$adj.P.Val < 0.1,
#'     down=function(df) df$logFC < 0 & df$adj.P.Val < 0.1
#' )
#'
#' processMicroarray(ma_data, treatment, control, classes, 'fdr')
#' processRNACounts(rna_data, treatment, control, classes, 'fdr')
#'
#' preprocessGeneExpressionData(dataList, classes, 'fdr')
#'
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes
#' @export
processMicroarray <- function(df, treatment, control, classes, adjust='none') {
   samples <- c(control, treatment)

   design <- matrix(0, length(samples), 2)
   colnames(design) <- c("control", "treatment")
   design[1:length(control), "control"] <- 1
   design[(length(control)+1):length(samples), "treatment"] <- 1
   contrastMatrix <- makeContrasts(DE=treatment-control, levels=design)

   fit <- lmFit(df, design)
   fit <- contrasts.fit(fit, contrastMatrix)
   fit <- eBayes(fit)

   .processFit(fit, classes, adjust)
}

#' @rdname preprocessDEGs
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom
#' @export
processRNACounts <- function(df, treatment, control, classes, adjust='none') {
   samples <- c(control, treatment)

   df <- df[, samples]

   design <- matrix(0, length(samples), 2)
   colnames(design) <- c("control", "treatment")
   design[1:length(control), "control"] <- 1
   design[(length(control)+1):length(samples), "treatment"] <- 1
   contrastMatrix <- makeContrasts(DE=treatment-control, levels=design)

   dge <- DGEList(df, remove.zeros = F)
   # TMM normalization
   dge <- calcNormFactors(dge)
   v <- voom(dge, design)

   fit <- lmFit(v, design)
   fit <- contrasts.fit(fit, contrastMatrix)
   fit <- eBayes(fit)

   .processFit(fit, classes, adjust)
}

#' @rdname preprocessDEGs
#' @export
preprocessGeneExpressionData <- function(dataList, classes, adjust='none') {
    common_genes <- Reduce(
        function(x, y) intersect(rownames(y$data), x),
        dataList[-1],
        rownames(dataList[[1]]$data)
    )
    if (length(common_genes) != nrow(dataList[[1]]$data)) {
        stop("Input data have different gene names")
    }
    geneIDs <- rownames(dataList[[1]]$data)
    sampleMatrix <- matrix(F, length(geneIDs), length(dataList),
                           dimnames=list(geneIDs, names(dataList)))

    result <- setNames(rep(list(sampleMatrix), length(classes)), names(classes))

    for (expName in names(dataList)) {
        item <- dataList[[expName]]
        if (item$type == 'MA') {
            df <- processMicroarray(item$data, item$treatment, item$control, classes, adjust)
        } else if (item$type == 'RNA') {
            df <- processRNACounts(item$data, item$treatment, item$control, classes, adjust)
        } else {
            stop("Wrong gene exprssion data type (must be 'MA' or 'RNA')")
        }

        for (name in names(classes)) {
            result[[name]][rownames(df), expName] <- df[, name]
        }
    }
    lapply(result, GeneClassificationMatrix)
}
