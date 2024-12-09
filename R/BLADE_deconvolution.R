#' Run BLADE deconvolution
#'
#' \code{BLADE_deconvolution.R} Runs BLADE to estimate cell fractions from bulk mRNA
#'
#'
#' @param signature SimpleList with mu and sigma, respectively mean gene
#' expression and mean-variance corrected variance per cell type
#' @param bulk mRNA to be deconvolved
#' @param genes subset of genes to be used for deconvolution
#' @param cores number of cores to use for paralellization
#' @param Alpha hyperparameter ....
#' @param Alpha0 hyperparameter ....
#' @param Kappa0 hyperparameter ....
#' @param SY hyperparameter ....
#' @param Nrep hyperparameter
#' @param Nrepfinal hyperparameter ....
#'
#' @return Statescope S4 object
#' @export
#'
#' @examples
#' predicted_fractions = BLADE_deconvolution(signature, bulk, selected_genes, 1L)
BLADE_deconvolution <- function(signature, bulk, genes, cores = 1L,
                      Alpha = 1L, Alpha0 = 1000L, Kappa0 = 1L, SY = 1L,
                      Nrep = 10L, Nrepfinal = 1000L) {
  ## save sample names for later use
  sample_names = colnames(bulk)

  ## init Mu en omega from list
  Mu = as.matrix(signature$mu[genes,])
  Omega = as.matrix(signature$omega[genes,])
  bulk = as.matrix(assay(bulk[genes, ],'normalized_counts'))

  result = Framework_Iterative(Mu, Omega,
                               bulk,
                               Alpha = Alpha, Alpha0= Alpha0,
                               Kappa0 = Kappa0, sY = SY, Nrep = Nrep,
                               Njob =cores, IterMax = Nrepfinal)

  ## make fractions Df
  fractions = DataFrame(t(result[[1]]$ExpF(result[[1]]$Beta)),
                        row.names = colnames(signature$mu))
  colnames(fractions) = sample_names

  ## Save S4 object with Statescope and fractions slots
  Statescope = new('Statescope', BLADE_output = result, fractions = fractions)

  return(Statescope)

}

