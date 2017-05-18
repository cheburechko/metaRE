#' @importFrom limma topTable
.processFit = function(fit, classes, adjust='none') {
   df <- topTable(fit, number=Inf, adjust.method=adjust)

   result <- data.frame(row.names=rownames(df))
   for (name in names(classes)) {
       result[ ,name] <- classes[[name]](df)
   }
   return(result)
}

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
