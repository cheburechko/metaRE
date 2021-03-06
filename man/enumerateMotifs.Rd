% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enumerateMotifs.R
\name{enumerateMotifs}
\alias{enumerateMotifs}
\alias{enumerateOligomers}
\alias{enumeratePatterns}
\title{Enumeration of various kinds of motifs.}
\usage{
enumerateOligomers(regulatoryRegions, k, rc = TRUE, output = c("genes",
  "counts", "positions", "composition"))

enumeratePatterns(regulatoryRegions, patterns, rc = TRUE,
  output = c("genes", "counts", "positions"))
}
\arguments{
\item{regulatoryRegions}{named charachter vector of nucleotide strings}

\item{k}{size of kmers}

\item{rc}{boolean, \code{TRUE} if motifs should be considered as equal to
their reverse complements (default \code{rc=TRUE})}

\item{output}{in which format the data should be returned (see Details)}

\item{patterns}{character vector of tested motifs}
}
\value{
Type of returned data structure depends on \code{output} parameter:
\describe{
    \item{genes}{named list of integer vectors, names are motifs, each
    vector describes a list of genes where the motif was found. Gene names
    are stored in 'genes' attribute of the list (default parameter). Only
    this format is recognised by \code{\link{calculateMassContingencyTablePvalues}}}
    \item{counts}{named integer vector, names are motifs, integers describe
    motif frequency in the given set}
    \item{positions}{named list of named lists, names are motifs, names in
    sublists are regulatory region names, each sublist describes exact
    positions of the last nucleotide of a motif in the corresponding region}
    \item{composition}{named list of character vectors, names are regulatory
    region names, each vector is an ordered sequence of oligomers that
    constitute this regulatory region}
}
}
\description{
Given a list of named regulatory regions, enumerate all possible
or only specific motifs and return data on their positions in these regions.
}
\details{
\code{enumerateOligomers} finds all possible oligomers of length
\code{k}.

\code{enumeratePatterns} looks only for specific motifs. Motifs can be
described in IUPAC nucleotide code with degenerate nucleotides.

Note that 'U' and '.' are not recognized by these functions.
}
\examples{
test_sequences <- c(
    gene1='aaaatgtcaaaa',
    gene2='ccccaaaagggg',
    gene3='ttttggggcccc'
)
k <- 4
enumerateOligomers(test_sequences, k, rc=FALSE)
enumerateOligomers(test_sequences, k, rc=TRUE)
enumerateOligomers(test_sequences, k, output='counts')
enumerateOligomers(test_sequences, k, output='positions')
enumerateOligomers(test_sequences, k, output='composition')
enumeratePatterns( test_sequences, c('AAAA', 'CCCC'), rc=FALSE, output='positions')
enumeratePatterns( test_sequences, c('AAAA', 'CCCC'), rc=TRUE, output='positions')
}
\seealso{
\code{\link{enumerateDyadsWithCore}}, \code{\link{enumerateRepeats}}
}
