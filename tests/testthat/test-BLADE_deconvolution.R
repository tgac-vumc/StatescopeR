library(scRNAseq)

test_that("BLADE deconvolution works properly with prior on simulation data", {
    ## Load scRNAseq
    scRNAseq <- scRNAseq::SegerstolpePancreasData()

    ## Preprocess scRNAseq
    scRNAseq$donor <- scRNAseq$individual
    scRNAseq$label <- scRNAseq$`cell type`

    ## Subset to 3 healthy and 3 type 2 diabetes samples
    scRNAseq = scRNAseq[,scRNAseq$donor %in% c('H2', 'H3',
                                              'T2D1', 'T2D2')]
    ## remove NA cells
    scRNAseq <- scRNAseq[, !is.na(scRNAseq$label)]

    ## remove cells with less than 100 in total cohort
    celltypes_to_remove <-
        names(table(scRNAseq$label)[(table(scRNAseq$label) < 100)])
    scRNAseq <- scRNAseq[, !scRNAseq$label %in% celltypes_to_remove]

    ##  pseudobulk
    load(system.file('extdata', 'example_pseudobulk.RData',
                     package = 'StatescopeR'))

    ##  Load signature
    load(system.file('extdata', 'example_signature.RData',
                     package = 'StatescopeR'))

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
