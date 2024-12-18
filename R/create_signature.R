#' Create Signature
#'
#' \code{create_signature} Creates signature from scRNAseq data for Deconvolution
#'
#'
#' @param data SingleCellExperiment object of which to make signature
#'
#' @return SimpleList DataFrames for Mu (mean per gene per cell type) and
#' Omega (variance corrected std.dev per gene per cell type)
#' @import matrixStats scran
#' @export
#'
#' @examples
#' signature = create_signature(data)
create_signature <- function(data) {
  ## init Mu, Omega & Var
  mu = DataFrame()
  omega = DataFrame()
  var = DataFrame()

  for (celltype in unique(data$label)){
    ## subset data on celltype
    temp_data = data[,data$label == celltype]
    ## Calculate Mu, Omega & Var
    ct_mu = rowMeans(as.array(logcounts(temp_data)))
    ct_omega = rowSds(as.array(logcounts(temp_data)))
    ct_var = rowVars(as.array(logcounts(temp_data)))

    ## Add ct Mu, Omega & var
    mu[celltype] = ct_mu
    omega[celltype] = ct_omega
    var[celltype] = ct_var

  }

  ## Correct Omega by mean-var trend
  new_omega <- omega[,0] ## make df with same # rows
  for (ct in colnames(omega)){
    trend <- fitTrendVar(mu[,ct], omega[,ct])$trend
    new_omega[,ct] <- trend(mu[,ct])

  }

  return(SimpleList(mu = mu, omega=new_omega))

}

