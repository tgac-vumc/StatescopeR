#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @import ComplexHeatmap
#' @importFrom stats cor setNames
#' @importFrom utils head
#'
NULL

#' Create a barplot of True vs estimated cfs
#'
#' Create a barplot of TRUE vs estimated cfs
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param true_fractions s4 Dataframe with true fractions to compare with
#'  estimated fractions
#' @return A barplot rendered to the active graphics device
#' @export fraction_eval
#' @examples
#' ## Load scRNAseq
#' scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' scRNAseq <- scRNAseq[1:100]
#' ## Preprocess scRNAseq
#' scRNAseq$donor <- scRNAseq$individual
#' scRNAseq$label <- scRNAseq$`cell type`
#'
#' ## remove NA cells
#' scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]
#'
#' ## remove duplicates gene names
#' scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove <-
#'     names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
#' scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]
#'
#' scRNAseq <- normalize_scRNAseq(scRNAseq)
#'
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk <- generate_pseudobulk(scRNAseq)
#'
#' pseudobulk <- normalize_bulkRNAseq(pseudobulk)
#'
#' ## Create signature from scRNAseq for deconvolution
#' signature <- create_signature(scRNAseq)
#'
#' ## Select genes optimized for deconvolution
#' selected_genes <- select_genes(scRNAseq)
#'
#' ## Optionally create prior expectation
#' prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
#' prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell
#'
#' ## Tranpose it to nSample x nCelltype
#' prior <- t(prior)
#'
#' ## Perform Deconvolution with BLADE, refine gene expression estimates
#' Statescope <- BLADE_deconvolution(
#'     signature, pseudobulk, selected_genes,
#'     prior, 1L
#' )
#'
#' ## ## Plot fraction correlation and RMSE per ct
#' fraction_eval(Statescope, gather_true_fractions(scRNAseq))
#'
fraction_eval <- function(Statescope, true_fractions) {
    ## measure correlation and RMSE per celltype
    eval_results <- setNames(
        data.frame(matrix(ncol = 3, nrow = 0)),
        c("celltype", "correlation", "RMSE"))
    for (ct in unique(rownames(true_fractions))) {
        cor_ct <- cor(
            as.matrix(true_fractions)[ct, ],
            as.matrix(fractions(Statescope))[ct, ])

        rmse_ct <- sqrt(mean((as.matrix(fractions(Statescope))[ct, ] -
            as.matrix(true_fractions)[ct, ])^2))

        ## add cor to cors
        eval_results[nrow(eval_results) + 1, ] <-
            data.frame(ct, cor_ct, rmse_ct)
    }

    ## plot correlation per celltype
    corplot <- ggplot(eval_results, aes(
        x = celltype, y = correlation,
        fill = celltype)) +
        geom_bar(stat = "identity", width = 0.95) +
        theme_bw() +
        labs(x = NULL) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

    ## plot correlation per celltype
    rmseplot <- ggplot(eval_results, aes(
        x = celltype, y = RMSE,
        fill = celltype)) +
        geom_bar(stat = "identity", width = 0.95) +
        theme_bw() +
        labs(x = NULL) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

    ## Plot fraction evaluation per ct
    plot_grid(corplot, rmseplot)
}

#' Create a heatmap of the estimated fractions
#'
#' Create a heatmap of the estimated fractions ...
#'
#' @param Statescope Statescope obj from BLADE_deconvolution
#' @param ... other parameters to pass to  [ComplexHeatmap::Heatmap()]
#' @return A heatmap rendered to the active graphics device
#' @export fraction_heatmap
#' @examples
#' #' ## Load scRNAseq
#' scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' scRNAseq <- scRNAseq[1:100]
#' ## Preprocess scRNAseq
#' scRNAseq$donor <- scRNAseq$individual
#' scRNAseq$label <- scRNAseq$`cell type`
#'
#' ## remove NA cells
#' scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]
#'
#' ## remove duplicates gene names
#' scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove <-
#'     names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
#' scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]
#'
#' scRNAseq <- normalize_scRNAseq(scRNAseq)
#'
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk <- generate_pseudobulk(scRNAseq)
#'
#' pseudobulk <- normalize_bulkRNAseq(pseudobulk)
#'
#' ## Create signature from scRNAseq for deconvolution
#' signature <- create_signature(scRNAseq)
#'
#' ## Select genes optimized for deconvolution
#' selected_genes <- select_genes(scRNAseq)
#'
#' ## Optionally create prior expectation
#' prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
#' prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell
#'
#' ## Tranpose it to nSample x nCelltype
#' prior <- t(prior)
#'
#' ## Perform Deconvolution with BLADE, refine gene expression estimates
#' Statescope <- BLADE_deconvolution(
#'     signature, pseudobulk, selected_genes,
#'     prior, 1L
#' )
#'
#' ## Plot fraction heatmap
#' fraction_heatmap(Statescope)
#'
fraction_heatmap <- function(Statescope, ...) {
    Heatmap(as.matrix(fractions(Statescope)),
        heatmap_legend_param = list(title = "")
    )
}


#' Create a barplot of top stateloadings
#'
#' Create a barplot of top stateloadings
#'
#' @param Statescope Statescope obj from StateDiscovery
#' @param top_n integer selecting how many genes to show per state
#' @return A barplot rendered to the active graphics device
#' @export barplot_stateloadings
#' @examples
#' #' ## Load scRNAseq
#' scRNAseq <- scRNAseq::SegerstolpePancreasData()
#'
#' ## subset to 100 genes for example
#' scRNAseq <- scRNAseq[1:100]
#' ## Preprocess scRNAseq
#' scRNAseq$donor <- scRNAseq$individual
#' scRNAseq$label <- scRNAseq$`cell type`
#'
#' ## remove NA cells
#' scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]
#'
#' ## remove duplicates gene names
#' scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
#'
#' ## remove cells with less than 100 in total cohort
#' celltypes_to_remove <-
#'     names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
#' scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]
#'
#' scRNAseq <- normalize_scRNAseq(scRNAseq)
#'
#' ## Create and normalized pseudobulk from scRNAseq
#' pseudobulk <- generate_pseudobulk(scRNAseq)
#'
#' pseudobulk <- normalize_bulkRNAseq(pseudobulk)
#'
#' ## Create signature from scRNAseq for deconvolution
#' signature <- create_signature(scRNAseq)
#'
#' ## Select genes optimized for deconvolution
#' selected_genes <- select_genes(scRNAseq)
#'
#' ## Optionally create prior expectation
#' prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
#' prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell
#'
#' ## Tranpose it to nSample x nCelltype
#' prior <- t(prior)
#'
#' ## Perform Deconvolution with BLADE, refine gene expression estimates
#' Statescope <- BLADE_deconvolution(
#'     signature, pseudobulk, selected_genes,
#'     prior, 1L
#' )
#'
#' ## Plot fraction heatmap
#' barplot_stateloadings(Statescope, top_n = 1)
#'
barplot_stateloadings <- function(Statescope, top_n = 1) {
    ## init df
    plot_df <- setNames(
        data.frame(matrix(ncol = 4, nrow = 0)),
        c("celltype", "state", "gene", "score")
    )
    ## Gather top genes and scores
    for (ct in names(stateloadings(Statescope))) {
        for (state in colnames(stateloadings(Statescope)[[ct]])) {
            indices <- head(sort(stateloadings(Statescope)[[ct]][, state],
                index.return = TRUE, decreasing = TRUE
            )$ix, top_n)
            genes <- rownames(stateloadings(Statescope)[[ct]][indices, ])
            scores <- stateloadings(Statescope)[[ct]][indices, state]

            ## add to plot df
            for (i in seq(top_n)) {
                plot_df[nrow(plot_df) + 1, ] <-
                    data.frame(ct, state, genes[i], scores[i])
            }
        }
    }

    ## Make barplot
    ggplot(plot_df, aes(x = celltype, y = score, fill = state)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_text(aes(label = gene),
            position = position_dodge(width = 0.9),
            vjust = 1.3
        )
}
