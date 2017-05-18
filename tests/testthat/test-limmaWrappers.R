library(testthat)

genes <- 1000
geneIDs <- paste0('gene', 1:genes)

experiment <- data.frame(
    control1=rep(4, genes),
    control2=rep(4, genes),
    control3=rep(4, genes),

    treatment1=c(rep(4, genes %/% 2), rep(1, genes %/% 4), rep(7, genes %/% 4)),
    treatment2=c(rep(4, genes %/% 2), rep(1, genes %/% 4), rep(7, genes %/% 4)),
    treatment3=c(rep(4, genes %/% 2), rep(1, genes %/% 4), rep(7, genes %/% 4)),

    row.names = geneIDs
)

rnaseqExperiment <- data.frame(
    control1=rep(16, genes),
    control2=rep(16, genes),
    control3=rep(16, genes),

    treatment1=c(rep(16, genes %/% 2), rep(2, genes %/% 4), rep(128, genes %/% 4)),
    treatment2=c(rep(16, genes %/% 2), rep(2, genes %/% 4), rep(128, genes %/% 4)),
    treatment3=c(rep(16, genes %/% 2), rep(2, genes %/% 4), rep(128, genes %/% 4)),

    row.names = geneIDs
)

control <- paste0('control', 1:3)
treatment <- paste0('treatment', 1:3)

classes <- list(
    up=function(df) df$logFC >= 1,
    up16=function(df) df$logFC >= 4,
    down=function(df) df$logFC <= -1,
    down16=function(df) df$logFC <= -4
)

answers <- data.frame(
    up=c(rep(F, genes %/% 4 * 3), rep(T, genes %/% 4)),
    up16=rep(F, genes),
    down=c(rep(F, genes %/% 2), rep(T, genes %/% 4), rep(F, genes %/% 4)),
    down16=rep(F, genes),
    row.names = geneIDs
)

test_that("microarrays are processed", {
    classification <- processMicroarray(experiment, treatment, control, classes)
    classification <- classification[rownames(answers), colnames(answers)]

    expect_equal(classification, answers)
})

test_that("rnaseqs are processed", {
    classification <- processRNACounts(rnaseqExperiment, treatment, control, classes)
    classification <- classification[geneIDs, ]

    expect_equal(classification, answers)
})


test_that("tasks are processed", {
    experiments <- list(
        exp1=list(data=experiment, treatment=treatment, control=control, type='MA'),
        exp2=list(data=rnaseqExperiment, treatment=treatment, control=control, type='RNA')
    )

    processedData <- preprocessGeneExpressionData(experiments, classes)

    expect_equal(names(processedData), names(classes))
    for (class in names(processedData)) {
        mat <- processedData[[class]][rownames(answers), ]
        expect_equal(colnames(mat), c('exp1', 'exp2'))
        for (column in colnames(mat)) {
            classification <- mat[, column]
            expect_equal(classification, setNames(answers[, class], rownames(answers)))
        }
    }



})
