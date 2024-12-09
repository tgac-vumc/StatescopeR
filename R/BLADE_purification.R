#' Run BLADE purification
#'
#' \code{BLADE_purification.R} Runs BLADE to estimate cell type-specific
#' gene expression profiles
#'
#' @param BLADE_result BLADE_result obj from BLADE_deconvolution
#' @param signature SimpleList with mu and sigma, respectively mean gene
#' expression and mean-variance corrected variance per cell type
#' @param bulk mRNA to be deconvolved
#' @param cores number of cores to use for paralellization
#'
#' @return updated BLADE_result S4 object with ct_specific_gep added
#' @export
#'
#' @examples
#' ct_specific_gep = BLADE_purification(BLADE_result, signature bulk, 1L)
BLADE_purification <- function(BLADE_result, signature, bulk, cores = 1L) {
  ## get final BLADE object from result ready for Python code
  BLADE_obj = list('final_obj' = BLADE_result@BLADE_output[[1]],
                   'outs' = BLADE_result@BLADE_output[[4]])

  ## init Mu en omega from list
  Mu = as.matrix(signature$mu)
  Omega = as.matrix(signature$omega)
  bulk = as.matrix(assay(pseudobulk,'normalized_counts'))

  result = Purify_AllGenes(BLADE_obj, Mu, Omega,
                               bulk, cores)

  ## update BLADE results
  BLADE_result@BLADE_output = result

  ## Gather ct specific gep
  ct_specific_gep = SimpleList()
  i = 0
  for (ct in colnames(signature$mu)){
    i = i+1
    temp_gep = t(result[[1]]$Nu[,,i])
    ## make absolute expression by fraction weighting
    temp_abs_gep = sweep(temp_gep, MARGIN =2,
                         as.matrix(BLADE_result@fractions[ct,]), `*`)
    ## Weight by Omega
    temp_omega_abs_gep = (temp_abs_gep - colMeans(temp_abs_gep)) *
                          result[[1]]$Omega[,i]

    ## Name samples and genes
    temp_omega_abs_gep_df = DataFrame(temp_omega_abs_gep)
    colnames(temp_omega_abs_gep_df) = colnames(bulk)
    rownames(temp_omega_abs_gep_df) = rownames(bulk)

    ## add to ct_specific_gp list
    ct_specific_gep[ct] =  SummarizedExperiment(assays = list(
      weighted_gep = temp_omega_abs_gep_df))
  }

  ## Add cell type specific gene expression to BLADE_result obj
  BLADE_result@ct_specific_gep = ct_specific_gep

  return(BLADE_result)

}

