library(scRNAseq)
library(StatescopeR)

test_that("BLADE deconvolution works properly with prior on simulation data", {
    ## Load SegerstolpePancreas data set
    scRNAseq <- SegerstolpePancreasData()
    scRNAseq$donor <- scRNAseq$individual
    scRNAseq$label <- scRNAseq$`cell type`

    ## remove cells with no cell type label
    scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]

    ## remove duplicate genes
    scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]

    ## remove very rare cell types (<100 cells in total data set)
    celltypes_to_remove <-
        names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
    scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]

    scRNAseq <- normalize_scRNAseq(scRNAseq)

    ## Create and normalized pseudobulk from scRNAseq
    pseudobulk <- generate_pseudobulk(scRNAseq)

    pseudobulk <- normalize_bulkRNAseq(pseudobulk)

    ## Create signature from scRNAseq for deconvolution
    signature <- create_signature(scRNAseq)

    ## Select genes optimized for deconvolution
    selected_genes <- select_genes(scRNAseq, 60L, n_hvg_genes = 200L)

    ## Optionally create prior expectation
    prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
    prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell

    ## Tranpose it to nSample x nCelltype
    prior <- t(prior)

    ## Perform Deconvolution with BLADE
    Statescope <- BLADE_deconvolution(
        signature, pseudobulk, selected_genes,
        prior, 2L
    )

    ## Compare true fractions with deconvolution results
    true_fractions = gather_true_fractions(scRNAseq)

    ## measure ct correlation with true fractions
    cors = list()
    for (ct in unique(rownames(true_fractions))){
        cor = cor(as.matrix(true_fractions)[ct,],
                  as.matrix(fractions(Statescope))[ct,])

        ## add cor to cors
        cors[ct]= cor

    }

    ## calculate median correlation with true fractions
    median_cor = median(unlist(cors))

    expect_gt(median_cor, 0.4)
})
