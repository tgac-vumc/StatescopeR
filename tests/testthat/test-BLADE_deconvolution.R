library(scRNAseq)

test_that("BLADE deconvolution works properly with prior on simulation data", {
    ## Load scRNAseq
    scRNAseq <- scRNAseq::SegerstolpePancreasData()

    ## remove duplicates gene names
    scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]

    ## Preprocess scRNAseq
    scRNAseq$donor <- scRNAseq$individual
    scRNAseq$label <- scRNAseq$`cell type`

    ## remove NA cells
    scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]

    ## remove cells with less than 100 in total cohort
    celltypes_to_remove <-
        names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
    scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]

    ## preprocessing
    scRNAseq <- normalize_scRNAseq(scRNAseq)

    ## Create and normalized pseudobulk from scRNAseq
    pseudobulk <- generate_pseudobulk(scRNAseq)

    pseudobulk <- normalize_bulkRNAseq(pseudobulk)

    ## Create signature from scRNAseq for deconvolution
    signature <- create_signature(scRNAseq, hvg_genes = TRUE,
                                  n_hvg_genes = 100L)

    ##  Load selected genes
    load(system.file('extdata', 'example_selected_genes.RData',
    package = 'StatescopeR'))

    ## Optionally create prior expectation
    prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
    prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell

    ## Tranpose it to nSample x nCelltype
    prior <- t(prior)

    ## Perform Deconvolution with BLADE, refine gene expression estimates
    Statescope <- BLADE_deconvolution(
        signature, pseudobulk, selected_genes,
        prior, 1L, Nrep = 1L ## Parallel causes workers to hang
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
