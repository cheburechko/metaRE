# metaRE
R package for motif discovery via meta-analysis of microarrays and RNA-Seq

## Installation

1. Install [Bioconductor](http://www.bioconductor.org/):
```R 
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
```

2. Install required packages:
```R
install.packages(c('Rcpp',  'BH', 'futile.logger', 'foreach'))
biocLite(c('limma', 'edgeR', 'GEOquery'))
```

3. Download metaRE:
```
git clone https://github.com/cheburechko/metaRE
```

4. Install metaRE:
From console
```
R CMD INSTALL [path/to/package]
```

or from R using [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

```R
install.packages('devtools')
devtools::install('[path/to/package]')
```
