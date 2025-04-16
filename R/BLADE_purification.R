#' Run BLADE purification
#'
#' \code{BLADE_purification.R} Runs BLADE to estimate cell type-specific
#' gene expression profiles
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param signature SimpleList with mu and sigma, respectively mean gene
#' expression and mean-variance corrected variance per cell type
#' @param bulk mRNA to be deconvolved
#' @param cores number of cores to use for paralellization
#'
#' @return updated Statescope S4 object with ct_specific_gep added
#' @import basilisk reticulate
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#'  ## Load data
#' data <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' data = data[1:100]
#' ## Preprocess data
#' data$donor = data$individual
#' data$label = data$`cell type`
#'
#' ## remove NA cells
#' data = data[,!is.na(data$label)]
#'
#' ## remove duplicates gene names
#' data = data[!duplicated(rownames(data)),]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove = names(table(data$label)[(table(data$label) <100)])
#' data = data[,!data$label %in% celltypes_to_remove]
#'
#' data = normalize_scRNAseq(data)
#'
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk = generate_pseudobulk(data)
#'
#' pseudobulk = normalize_bulkRNAseq(pseudobulk)
#'
#' ## Measure true cell fractions from pseudobulk
#' true_fractions = gather_true_fractions(data)
#'
#' ## Create signature from scRNAseq for deconvolution
#' signature = create_signature(data)
#'
#' ## Select genes optimized for deconvolution
#' selected_genes = select_genes(data)
#'
#' ## Perform Deconvolution with BLADE
#' Statescope = BLADE_deconvolution(signature, pseudobulk, selected_genes, 1L)
#' Statescope = BLADE_purification(Statescope, signature, pseudobulk, 1L)
#' ct_specific_getp = Statescope@ct_specific_gep
BLADE_purification <- function(Statescope, signature, bulk, cores = 1L) {
  ## get final BLADE object from result ready for Python code
  BLADE_obj = list('final_obj' = Statescope@BLADE_output[[1]],
                   'outs' = Statescope@BLADE_output[[3]])

  ## init Mu en omega from list
  Mu = as.matrix(signature$mu)
  Omega = as.matrix(signature$omega)
  bulk = as.matrix(assay(bulk,'normalized_counts'))

  ## start basilisk
  setBasiliskShared(TRUE)
  proc <- basiliskStart(deconvolution)

  ## Refine gene expression estimation using Basilisk
  Statescope <- basiliskRun(proc, fun = function(BLADE_obj, Mu, Omega,
                                             bulk, cores){
      ## import OncoBLADE
      reticulate::source_python(system.file('python/OncoBLADE.py',
                                            package ='StatescopeR'))

      ## Run refinement
      result = Purify_AllGenes(BLADE_obj, Mu, Omega,
                               bulk, cores)

      ## source StatescopeR package for class and DataFrame for fractions
      library(StatescopeR)
      #library(S4Vectors)

      ## update BLADE results
      Statescope@BLADE_output = result

      ## Gather ct specific gep
      ## TODO:
      ## -Find solution for this with SimpleList or SummarizedExperiment
      ct_specific_gep = list()
      i = 0
      for (ct in colnames(signature$mu)){
          i = i+1
          temp_gep = t(Statescope@BLADE_output[[1]]$Nu[,,i])
          ## make absolute expression by fraction weighting
          temp_abs_gep = sweep(temp_gep, MARGIN =2,
                               as.matrix(Statescope@fractions[ct,]), `*`)
          ## Weight by Omega
          temp_omega_abs_gep = (temp_abs_gep - colMeans(temp_abs_gep)) *
              Statescope@BLADE_output[[1]]$Omega[,i]

          ## Name samples and genes
          temp_omega_abs_gep_df = DataFrame(temp_omega_abs_gep)
          colnames(temp_omega_abs_gep_df) = colnames(bulk)
          rownames(temp_omega_abs_gep_df) = rownames(bulk)

          ## add to ct_specific_gp list
          ct_specific_gep[ct] =  SummarizedExperiment(assays = SimpleList(
              weighted_gep = temp_omega_abs_gep_df))
      }

      ## Add cell type specific gene expression to Statescope obj
      Statescope@ct_specific_gep = ct_specific_gep
      Statescope

  }, BLADE_obj=BLADE_obj , Mu=Mu, Omega=Omega, bulk=bulk, cores=cores)

  ## stop basilisk
  basiliskStop(proc)

  return(Statescope)

}

