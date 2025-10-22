## usethis namespace: start
## usethis namespace: end
NULL

#' The StatescopeR package
#'
#' framework for discovery of cell states from bulk mRNA profiles
#'
#' @rdname StatescopeR-package
#' @name StatescopeR-package
#' @keywords internal
#' @aliases StatescopeR-package StatescopeR
#'
#' @section Details:
#'  StatescopeR starts from lognormalized cp10k single cell data and a bulk RNA
#'  dataset to be used for deconvolution and cell state analysis. From this
#'  point StatescopeR has some functions to create a signature and select genes
#'  for deconvolution.
#'  \code{\link{create_signature}} creates a signature of lognormalized cp10k
#'      single cell data for deconvolution and cell state analysis.
#'  \code{\link{select_genes}} selects genes which best discriminate cell types
#'      for use in deconvolution.
#'  \code{\link{BLADE_deconvolution}} estimates cell fractions from bulk mRNA.
#'  \code{\link{Refinement}} refines cell type-specific gene expression profile
#'      estimates to better capture inter-sample variability.
#'  \code{\link{StateDiscovery}} discovery cell states from inferred cell
#'  type-specific gene expression profiles.
#'
#'  The following functions provide some evaluation and plotting functions:
#'  \code{\link{gather_true_fractions}} gets true cell fractions from
#'      scRNAseq data.
#'  \code{\link{fraction_eval}} plots the correlation of estimated vs true
#'      cell fractionsover all samples.
#'  \code{\link{fraction_heatmap}} plots a heatmap of estimated
#'      cell fractions.
#'  \code{\link{barplot_stateloadings}} plots a barplot showing the top n
#'      most important stateloadings.
#'

"_PACKAGE"
