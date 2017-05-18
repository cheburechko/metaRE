#' @name enumerateMotifs
#' @title Enumeration of various kinds of motifs.
#' @param regulatoryRegions named charachter vector of nucleotide strings
#' @param k size of kmers
#' @param rc boolean, \code{TRUE} if motifs should be considered as equal to
#' their reverse complements (default \code{rc=TRUE})
#' @param output in which format the data should be returned (see Details)
#' @param patterns character vector of tested motifs
#' @description Given a list of named regulatory regions, enumerate all possible
#' or only specific motifs and return data on their positions in these regions.
#' @details \code{enumerateOligomers} finds all possible oligomers of length
#' \code{k}.
#'
#' \code{enumeratePatterns} looks only for specific motifs. Motifs can be
#' described in IUPAC nucleotide code with degenerate nucleotides.
#'
#' Note that 'U' and '.' are not recognized by these functions.
#' @seealso \code{\link{enumerateDyadsWithCore}}
#' @return Type of returned data structure depends on \code{output} parameter:
#' \describe{
#'     \item{genes}{named list of integer vectors, names are motifs, each
#'     vector describes a list of genes where the motif was found. Gene names
#'     are stored in 'genes' attribute of the list (default parameter). Only
#'     this format is recognised by \code{\link{calculateMassContingencyTablePvalues}}}
#'     \item{counts}{named integer vector, names are motifs, integers describe
#'     motif frequency in the given set}
#'     \item{positions}{named list of named lists, names are motifs, names in
#'     sublists are regulatory region names, each sublist describes exact
#'     positions of the last nucleotide of a motif in the corresponding region}
#'     \item{composition}{named list of character vectors, names are regulatory
#'     region names, each vector is an ordered sequence of oligomers that
#'     constitute this regulatory region}
#' }
#' @examples
#' test_sequences <- c(
#'     gene1='aaaatgtcaaaa',
#'     gene2='ccccaaaagggg',
#'     gene3='ttttggggcccc'
#' )
#' k <- 4
#' enumerateOligomers(test_sequences, k, rc=FALSE)
#' enumerateOligomers(test_sequences, k, rc=TRUE)
#' enumerateOligomers(test_sequences, k, output='counts')
#' enumerateOligomers(test_sequences, k, output='positions')
#' enumerateOligomers(test_sequences, k, output='composition')
#' enumeratePatterns( test_sequences, c('AAAA', 'CCCC'), rc=FALSE, output='positions')
#' enumeratePatterns( test_sequences, c('AAAA', 'CCCC'), rc=TRUE, output='positions')
NULL

.enumerateMotifs <- function(parameters) {
    enumerateMotifsCpp(parameters, GeneClassificationSparse, futile.logger::flog.debug)
}

#' @rdname enumerateMotifs
#' @export
enumerateOligomers <- function(regulatoryRegions, k, rc=TRUE,
                               output=c('genes', 'counts', 'positions', 'composition')) {
    .enumerateMotifs(list(
        regulatoryRegions=regulatoryRegions,
        counter=list(mode='simple', k=k, rc=rc), data=match.arg(output)
    ))
}

#' @name enumerateDyadsWithCore
#' @title Enumerate Dyads With Predefined Core
#' @param regulatoryRegions named charachter vector of nucleotide strings
#' @param k size of kmers
#' @param core character vector of possible core motifs in a dyad
#' @param minSpacer minimal distance in base pairs between core and a second
#' motif
#' @param maxSpacer maximal distance in base pairs between core and a second
#' motif
#' @param rc boolean, \code{TRUE} if motifs should be considered as equal to
#' their reverse complements (default \code{rc=TRUE})
#' @param output in which format the data should be returned
#' (see \code{\link{enumerateMotifs}})
#' @param fuzzySpacer if \code{TRUE}, dyads with the same core and partner
#' oligomer but different spacers will be counted as the same dyad.
#' @param fuzzyOrder if \code{TRUE}, dyads with the same core, partner oligomer
#' and spacer but different order (core is left/right part of dyad) will be
#' counted as the same dyad.
#' @param fuzzyOrientation if \code{TRUE} and \code{rc=TRUE} then dyads with
#' partner and its reverse complement will be counted as the same dyad.
#' @description Given a list of named regulatory regions, enumerate all possible
#' spaced dyads with a gicen core located within the defined spacer range
#' and return data on their positions in these regions.
#' @seealso \code{\link{enumerateMotifs}}
#' @examples
#' test_sequences <- c(
#'     gene1='ccccggggtgtcaaaccccc'
#' )
#' k <- 4
#' interesting_elements <- c('TGTC_4_CCCC', 'TGTC_3_CCCC', 'CCCC_4_TGTC')
#' result <- enumerateDyadsWithCore(test_sequences, k, 'TGTC', 0, 4,
#'                                        output='positions')
#' result[interesting_elements]
#' result <- enumerateDyadsWithCore(test_sequences, k, 'TGTC', 0, 4,
#'                                        fuzzyOrder=TRUE, output='positions')
#' result[interesting_elements]
#' interesting_elements <- c('TGTC_0..4_CCCC', 'CCCC_0..4_TGTC')
#' result <- enumerateDyadsWithCore(test_sequences, k, 'TGTC', 0, 4,
#'                                        fuzzySpacer=TRUE, output='positions')
#' result[interesting_elements]
#' result <- enumerateDyadsWithCore(test_sequences, k, 'TGTC', 0, 4,
#'                                        fuzzySpacer=TRUE, fuzzyOrder=TRUE,
#'                                        output='positions')
#' result[interesting_elements]
#' result <- enumerateDyadsWithCore(test_sequences, k, 'TGTC', 0, 4,
#'                                        fuzzySpacer=TRUE, fuzzyOrder=TRUE,
#'                                        fuzzyOrientation=TRUE,
#'                                        output='positions')
#' result['TGTC_0..4_GGGG']
#' @export
enumerateDyadsWithCore <- function(regulatoryRegions, k, core,
    minSpacer, maxSpacer, rc=TRUE, output=c('genes', 'counts', 'positions'),
    fuzzySpacer=FALSE, fuzzyOrder=FALSE, fuzzyOrientation=FALSE) {
    .enumerateMotifs(list(
        regulatoryRegions=regulatoryRegions, counter=list(
            mode='spaced_dyad', patterns=core, k=k, rc=rc, maxSpacer=maxSpacer,
            minSpacer=minSpacer, fuzzySpacer=fuzzySpacer, fuzzyOrder=fuzzyOrder,
            fuzzyOrientation=fuzzyOrientation
        ), data=match.arg(output)
    ))
}

#' @rdname enumerateMotifs
#' @export
enumeratePatterns <- function(regulatoryRegions, patterns, rc=TRUE,
                              output=c('genes', 'counts', 'positions')) {
    .enumerateMotifs(list(
        regulatoryRegions=regulatoryRegions,
        counter=list(mode='specific_single', patterns=patterns, rc=rc),
        data=match.arg(output)
    ))
}

#' @export
enumerateRepeats <- function(regulatoryRegions, k, minSpacer, maxSpacer,
                             rc=TRUE, output=c('genes', 'counts', 'positions'))
{
    .enumerateMotifs(list(
        regulatoryRegions=regulatoryRegions, counter=list(
            mode='repeat', k=k, rc=rc, maxSpacer=maxSpacer,
            minSpacer=minSpacer), data=match.arg(output)
    ))
}
