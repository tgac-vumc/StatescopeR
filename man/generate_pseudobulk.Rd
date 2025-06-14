% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_pseudobulk.R
\name{generate_pseudobulk}
\alias{generate_pseudobulk}
\title{Generate pseudobulk from single cell RNAseq}
\usage{
generate_pseudobulk(data)
}
\arguments{
\item{data}{SingleCellExperiment object of which to generate pseudobulk}
}
\value{
SummarizedExperiment with pseudobulk from the scRNAseq
}
\description{
\code{generate_pseudobulk} Creates pseudobulk from raw counts of
SingleCellExperiment object
}
\examples{
#' ## Load data
data <- scRNAseq::SegerstolpePancreasData()

## subset to 100 genes for example
data <- data[1:100]
## Preprocess data
data$donor <- data$individual
data$label <- data$`cell type`
pseudobulk <- generate_pseudobulk(data)
}
