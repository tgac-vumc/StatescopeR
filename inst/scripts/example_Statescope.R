# This script was used to create the `Statescope_Deconvolved.RData`,
# `Statescope_Refined.RData` and `Statescope_Discovered.RData` files.
# This file contains an example Statescope object for use in examples and tests.
# For this purpose the SegerstolpePancreasData from the scRNAseq package was used.
# This dataset was subset to the first 100 genes of all cells with celltype labels.
# Rare celltypes were excluded, leaving 6 celltypes for analysis.
# After this selection procedure, standard preprocessing was done,
# before running the Statescope framework and saving all steps.
# Key package versions:
#
# scRNAseq   v2.23.0

library(StatescopeR)
library(scRNAseq)

## Load scRNAseq
scRNAseq <- scRNAseq::SegerstolpePancreasData()

## remove duplicates gene names
scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]

## subset to 100 genes for example
scRNAseq <- scRNAseq[1:100]
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
signature <- create_signature(scRNAseq)

## Select genes optimized for deconvolution (small number of genes for speed)
selected_genes <- select_genes(scRNAseq, 100L)

## Optionally create prior expectation
prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell

## Tranpose it to nSample x nCelltype
prior <- t(prior)

## Perform Deconvolution with BLADE, refine gene expression estimates
Statescope <- BLADE_deconvolution(
    signature, pseudobulk, selected_genes,
    prior, 2L
)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Deconvolved.RData')

Statescope <- Refinement(Statescope, signature, pseudobulk, 2L)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Refined.RData')

## Discover states
Statescope <- StateDiscovery(Statescope, Ncores = 2L, max_clusters = 4L)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Discovered.RData')

