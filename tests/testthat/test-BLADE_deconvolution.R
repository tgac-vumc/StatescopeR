library(scRNAseq)
library(StatescopeR)

test_that("BLADE deconvolution works properly  with prior on simulation data", {
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

    ## Normalize (cp10k) and logtransform scRNAseq
    scRNAseq <- normalize_scRNAseq(scRNAseq)

    ## Create scRNAseq reference/signature
    signature <- create_signature(scRNAseq, hvg_genes = TRUE)

    ## select subset of genes for deconvolution
    selected_genes <- select_genes(scRNAseq, 200L) # 200 genes to make it quick

    ## Create pseudobulk and also lognormalize
    pseudobulk <- generate_pseudobulk(scRNAseq)
    pseudobulk <- normalize_bulkRNAseq(pseudobulk)

    ## (optional) Create prior expectation
    prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
    prior[rownames(prior) != "ductal cell", ] <- NA #Keep only ductal cells as prior
    prior <- t(prior) # Tranpose it to nSample x nCelltype

    ## Run Deconvolution module
    Statescope <- BLADE_deconvolution(signature, pseudobulk, selected_genes, prior,
                                      cores = 1L)

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

    expect_gt(median_cor, 0.6)
})
