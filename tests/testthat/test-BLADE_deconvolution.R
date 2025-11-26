library(scRNAseq)
library(scuttle)

test_that("BLADE deconvolution works properly with prior on simulation data", {
    ## Load SegerstolpePancreas data set
    scRNAseq <- SegerstolpePancreasData()

    ## remove duplicate genes
    scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]
    ## Subset to 2 healthy and 3 type 2 diabetes samples
    scRNAseq = scRNAseq[,scRNAseq$individual %in% c('H2', 'H3',
                                                    'T2D1', 'T2D2')]
    ## remove cells with no cell type label
    scRNAseq <- scRNAseq[, !is.na(scRNAseq$`cell type`)]
    ## remove very rare cell types (<120 cells in total data set)
    celltypes_to_remove <-names(table(scRNAseq$`cell type`)
             [(table(scRNAseq$`cell type`) < 120)])
    scRNAseq <- scRNAseq[, !scRNAseq$`cell type` %in% celltypes_to_remove]

    ## Subset to first 3k genes
    scRNAseq = scRNAseq[1:3000,]

    ## Create pseudobulk and normalize to cp10k
    pseudobulk <- aggregateAcrossCells(scRNAseq, ids = scRNAseq$individual)
    normcounts(pseudobulk) <- calculateCPM(pseudobulk)/100
    pseudobulk = as(pseudobulk, "SummarizedExperiment")
    rownames(pseudobulk) = rownames(scRNAseq)

    ##  Load signature
    load(system.file("extdata", "example_signature.RData",
        package = "StatescopeR"
    ))

    ##  Load selected genes
    load(system.file("extdata", "example_selected_genes.RData",
        package = "StatescopeR"
    ))

    ##  Load prior
    load(system.file("extdata", "example_prior.RData",
        package = "StatescopeR"
    ))

    ## Perform Deconvolution with BLADE, refine gene expression estimates
    Statescope <- BLADE_deconvolution(
        signature, pseudobulk, selected_genes,
        prior, 1L,
        Nrep = 1L ## Parallel causes workers to hang
    )

    ##  Load true fractions
    load(system.file("extdata", "example_true_fractions.RData",
        package = "StatescopeR"
    ))

    ## measure ct correlation with true fractions
    cors <- list()
    for (ct in unique(rownames(true_fractions))) {
        cor <- cor(
            as.matrix(true_fractions)[ct, ],
            as.matrix(metadata(Statescope)$fractions)
            [ct, names(true_fractions)]
        )

        ## add cor to cors
        cors[ct] <- cor
    }

    ## calculate median correlation with true fractions
    median_cor <- median(unlist(cors))

    expect_gt(median_cor, 0.4)
})
