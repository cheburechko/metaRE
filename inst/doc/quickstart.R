## ---- eval=FALSE---------------------------------------------------------
#  library(metaRE)
#  MA_AT_auxin

## ---- eval=FALSE---------------------------------------------------------
#  experiments <- list()
#  
#  for (ma in MA_AT_auxin) {
#      # Download the data from GEO, convert to log2 if necessary.
#      exprs <- prepareGEO(control=ma$control, treatment=ma$treatment,
#                           isLog2=ma$log2)
#  
#      # Store in experiments for further processing
#      # Here 'type' can be either 'MA' or 'RNA', which tells if dataset is
#      # microarray expression levels or RNA-Seq RNA counts.
#      experiments[[ma$name]] <- list(
#          control=ma$control, treatment=ma$treatment, type="MA", data=exprs
#      )
#  }

## ---- eval=FALSE---------------------------------------------------------
#  # We will analyze up- and down-regulation separately. We will define 2
#  # functions that will select up- and down- regulated genes in 'classes' varaible.
#  logFCThreshold <- log2(1.5)
#  pValThreshold <- 0.05
#  classes <- list(
#      up=function(df) df$logFC > logFCThreshold & df$adj.P.Val < pValThreshold,
#      down=function(df) df$logFC < -logFCThreshold & df$adj.P.Val < pValThreshold
#  )
#  
#  # This will create a list of GeneClassificationMatix objects, one matrix per DEG
#  # class. GeneClassificationMatix is a logical matrix, where each row corresponds
#  # to a gene, each column corresponds to a dataset, and values tell if a gene
#  # belongs to this DEG class in this dataset.
#  # The last parameter is mulitple testing correction method as accepted by p.adjust
#  degs <- preprocessGeneExpressionData(experiments, classes, 'fdr')
#  head(degs$up)
#  head(degs$down)

## ---- eval=FALSE---------------------------------------------------------
#  ma <- experiments[[1]]
#  # The syntax is exactly the same for processRNACounts
#  singleExperiment <- processMicroarray(ma$data, ma$treatment, ma$control, classes, 'fdr')
#  
#  head(singleExperiment)

## ---- eval=FALSE---------------------------------------------------------
#  # Let 'up', 'down' be your custom DEG matrix
#  genes <- 1000
#  exps <- 6
#  total <- genes*exps
#  
#  up <- matrix(runif(total) < 0.05, genes, exps)
#  down <- matrix(runif(total) < 0.05, genes, exps)
#  custom_degs <- list(
#      up=GeneClassificationMatrix(up),
#      down=GeneClassificationMatrix(down)
#  )

## ---- eval=FALSE---------------------------------------------------------
#  promoters <- Biostrings::readDNAStringSet('/path/to/file.fasta')
#  promoters <- setNames(as.character(promoters), names(promoters))

## ---- eval=FALSE---------------------------------------------------------
#  library(biomaRt)
#  # Connect to Arabidopsis thaliana mart.
#  mart <- useMart('plants_mart', host="plants.ensembl.org", dataset='athaliana_eg_gene')
#  # Get gene promoters by ATH1 identifiers from our DEG matrix.
#  promoters <- getBM(
#      attributes=c("affy_ath1_121501", "gene_flank"),
#      filters=c("affy_ath1_121501","upstream_flank"),
#      values=list(rownames(degs$up), 1500),
#      mart=mart,
#      checkFilters=FALSE,
#      bmHeader=TRUE
#  )
#  # Convert data.frame to named character vector
#  promoters <- setNames(promoters$`Flank (Gene)`, promoters$`AFFY ATH1 121501 probe`)

## ---- eval=FALSE---------------------------------------------------------
#  # This function finds all possible oligomers of the fixed length in the promoter
#  # vector and returns a GeneClassificationSparse object.
#  # By default rc=TRUE which means that oligomers are considered equal to their
#  # reverse complements (e.g. AACCGG == CCGGTT)
#  regElements <- enumerateOligomers(promoters, 6)
#  
#  # The result is a named list of integer vectors. Names are cis-regulatory
#  # elements, vectors are indices of genes in which this elements are present.
#  # Gene names are stored in 'geneNames' attribute of the result, which can be
#  # accessed by 'geneNames(x)' funciton.
#  head(regElements)

## ---- eval=FALSE---------------------------------------------------------
#  data <- list(
#      elem1=c(1, 5, 10),
#      elem2=c(3, 6, 9),
#      elem3=2
#  )
#  genes <- paste("gene", 1:10)
#  regElements <- GeneClassificationSparse(data, genes)

## ---- eval=FALSE---------------------------------------------------------
#  pvalues <- list()
#  
#  for (class in names(degs)) {
#      # By default this looks for enrichment.
#      # Use "alternative='less'" for depletion or "alternative='two.sided'" for both.
#      pvalues[[class]] <- calculateMassContingencyTablePvalues(regElements, degs[[class]])
#  
#      # The result is a float matrix of p-values, where rows correspond to
#      # cis-regulatory elements, columns correspond to datasets, values are raw
#      # p-values.
#      print(class)
#      head(pvalues[[class]])
#  }
#  

## ---- eval=FALSE---------------------------------------------------------
#  metaPvalues <- list()
#  
#  
#  for (class in names(pvalues)) {
#      # 'adjust' is multiple testing correction method
#      # 'threshold' is a cutoff value for adjusted meta p-value
#      metaPvalues[[class]] <- calcMetaAssociation(
#          pvalues[[class]], adjust='bonferroni', threshold=0.05
#      )
#  
#      # The result is a data.frame with Meta.P.Value and Adj.Meta.P.Value columns.
#      # All elements with Adj.Meta.P.Value above cutoff threshold are removed from
#      # results
#      print(class)
#      head(metaPvalues[[class]])
#  }

## ---- eval=FALSE---------------------------------------------------------
#  metaPvalues <- list()
#  
#  for (class in names(degs)) {
#      # This performs both 'calculateMassContingencyTablePvalues' and
#      # 'calcMetaAssociation' in order.
#      metaPvalues[[class]] <- testRegulationHypotheses(
#          regElements, degs[[class]], adjust='bonferroni', threshold=0.05
#      )
#  }

## ---- eval=FALSE---------------------------------------------------------
#  # Enable doParallel backend for parallel processing of permutation test.
#  library(doParallel)
#  registerDoParallel()
#  
#  results <- list()
#  for (class in names(metaPvalues)) {
#      # Analyze only elements that pass the meta p-value threshold.
#      sigRegElements <- GeneClassificationSparse(
#          regElements[rownames(metaPvalues[[class]])],
#          geneNames(regElements)
#      )
#  
#      # Calculate the number of permutations required to obtain reasonable
#      # p-values to pass the Bonferroni adjusted cutoff threshold
#      threshold <- 0.05
#      n <- integer(length(sigRegElements)*50/threshold)
#  
#      # Permuation test every 'perRun' iterations stores the preliminary results
#      # in 'outfile' and removes elements that in future will not be able to
#      # pass the 'pvaluePreFilter' threshold.
#      outfile <- tempfile()
#      results[[class]] <- metaRE::permutationTest(
#          sigRegElements, degs[[degClass]], n, outfile=outfile,
#          pvaluePreFilter=threshold/length(sigRegElements), perRun=5000
#      )
#  
#      print(class)
#      head(results[[class]])
#  
#      # After the permutation test you can select the elements that pass the
#      # threshold and consider them as significantly associated with stimulus
#      # response.
#      print(results[[class]]$Hypothesis[results[[class]]$Permutation.P.Value < threshold/length(sigRegElements)])
#  }
#  

