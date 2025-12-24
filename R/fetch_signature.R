#' Fetch scRNAseq Signature
#'
#' \code{fetch_signature} Fetches signature from StatescopeData repository
#'
#'
#' @param tumor_type string choosing the tumor type from available signatures
#'  at: https://github.com/tgac-vumc/StatescopeData
#' @param n_celltypes integer choosing the number of cell types of the signature
#'  from: https://github.com/tgac-vumc/StatescopeData
#'
#' @return SimpleList with DataFrames for Mu (mean per gene per cell type) and
#' Omega (variance corrected std.dev per gene per cell type) and vector of
#' selected genes for deconvolution chosen by AutoGeneS
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' signature <- fetch_signature('PDAC', 7)
#' selected_genes <- signature$selected_genes
fetch_signature <- function(tumor_type, n_celltypes) {
    ## Check input formatting
    if (!is(tumor_type, "character")) {
        stop("tumor_type is not a string")
        } else if(!is(n_celltypes, "numeric")) {
            stop("n_celltypes is not numeric")}

    # Construct file name and URL
    file_name <- paste0(tumor_type, "_Signature_", n_celltypes, "celltypes.txt")
    file_url <- paste0(
        "https://raw.githubusercontent.com/tgac-vumc/StatescopeData/master/",
        tumor_type, "/", file_name)

    # Try to read the file; handle errors
    tryCatch({con <- url(file_url, open = "r")
            on.exit(close(con), add = TRUE)  # Ensure the connection closes
            df <- read.table(
                con,
                sep = "\t",
                header = TRUE,
                row.names = "Gene",
                check.names = FALSE)},
        error = function(e) {
            stop("TumorType not specified or invalid. Check
                https://github.com/tgac-vumc/StatescopeData for valid options.",
                call. = FALSE
                )
        })

    # Adjust to Bioconductor file formats
    mu_cols <- grep('scExp', colnames(df))
    omega_cols <- grep('scVar', colnames(df))

    mu <- S4Vectors::DataFrame(df[mu_cols])
    omega <- S4Vectors::DataFrame(df[omega_cols])
    selected_genes <- rownames(df)[as.logical(df$IsMarker)]

    ## Adjust colnames of mu and omega
    colnames(mu) <- sub("^.*_", "", colnames(mu))
    colnames(omega) <- sub("^.*_", "", colnames(omega))

    return(SimpleList(mu = mu, omega = omega, selected_genes = selected_genes))
}
