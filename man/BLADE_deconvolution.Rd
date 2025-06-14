% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLADE_deconvolution.R
\name{BLADE_deconvolution}
\alias{BLADE_deconvolution}
\title{Run BLADE deconvolution}
\usage{
BLADE_deconvolution(
  signature,
  bulk,
  genes,
  prior = NULL,
  cores = 1L,
  Alpha = 1L,
  Alpha0 = 1000L,
  Kappa0 = 1L,
  sY = 1L,
  Nrep = 10L,
  Nrepfinal = 1000L
)
}
\arguments{
\item{signature}{SimpleList with mu and sigma, both nGene x nCelltype
dataframes, respectively mean gene expression and mean-variance corrected
variance per cell type}

\item{bulk}{nSample x nGene mRNA to be deconvolved}

\item{genes}{subset of genes to be used for deconvolution}

\item{prior}{(optional) nSample x nCelltype matrix with prior fraction
expectations}

\item{cores}{number of cores to use for paralellization}

\item{Alpha}{BLADE Hyperparameter}

\item{Alpha0}{BLADE Hyperparameter}

\item{Kappa0}{BLADE Hyperparameter}

\item{sY}{BLADE Hyperparameter}

\item{Nrep}{Number of BLADE initializations}

\item{Nrepfinal}{Number of maximum optimization iterations}
}
\value{
Statescope S4 object
}
\description{
\code{BLADE_deconvolution.R} Runs BLADE to estimate cell fractions from bulk
mRNA
}
\examples{
## Load scRNAseq
scRNAseq <- scRNAseq::SegerstolpePancreasData()

## subset to 100 genes for example
scRNAseq <- scRNAseq[1:100]
## Preprocess scRNAseq
scRNAseq$donor <- scRNAseq$individual
scRNAseq$label <- scRNAseq$`cell type`

## remove NA cells
scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]

## remove duplicates gene names
scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]

## remove cells with less than 100 in total cohort
celltypes_to_remove <-
    names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
scRNAseq <- scRNAseq[, !scRNAseq$label \%in\% celltypes_to_remove]

scRNAseq <- normalize_scRNAseq(scRNAseq)

## Create and normalized pseudobulk from scRNAseq
pseudobulk <- generate_pseudobulk(scRNAseq)

pseudobulk <- normalize_bulkRNAseq(pseudobulk)

## Create signature from scRNAseq for deconvolution
signature <- create_signature(scRNAseq)

## Select genes optimized for deconvolution
selected_genes <- select_genes(scRNAseq)

## Optionally create prior expectation
prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell

## Tranpose it to nSample x nCelltype
prior <- t(prior)

## Perform Deconvolution with BLADE
Statescope <- BLADE_deconvolution(
    signature, pseudobulk, selected_genes,
    prior, 1L
)
fractions(Statescope)

}
