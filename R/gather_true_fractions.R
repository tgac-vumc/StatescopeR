#' Gather true fractions
#'
#' \code{gather_true_fractions} Gathers true fractions of all cell types per sample from scRNAseq data
#'
#'
#' @param data SingleCellExperiment object of which to gather fractions per sample
#'
#' @return DataFrame with fractions of all cell types per sample
#' @export
#'
#' @examples
#' true_fractions = gather_true_fractions(data)
gather_true_fractions <- function(data) {
  ## Init list to save fractions
  true_fractions = list()

  ## loop over samples
  for (sample in unique(colData(data)$donor)) {
    temp_data = data[,colData(data)$donor == sample]
    temp_true_fractions = DataFrame(table(temp_data$label)/ncol(temp_data))

    ## add True_fractions of sample to list
    true_fractions[sample] = list(temp_true_fractions$Freq)
  }

  ## Add all samples to one DataFrame
  true_fractions = DataFrame(true_fractions,
                             row.names = temp_true_fractions$Var1)

  ## return DF
  return(true_fractions)

}

