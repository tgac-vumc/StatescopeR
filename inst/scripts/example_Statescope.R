# This script was used to create the `example_scRNAseq.RData`, `example_signature.RData`,
# `example_selected_genes.RData`, `example_Statescope_Deconvolved.RData`,
# `example_Statescope_Refined.RData` and `example_Statescope_Discovered.RData` files.
# This file contains an example Statescope pipeline of which the intermediate results
# are saved for use in examples and tests.
# For this purpose the SegerstolpePancreasData from the scRNAseq package was used.
# Rare celltypes were excluded, leaving the 5 most common celltypes for analysis.
# After this exclusion, standard preprocessing was done, after which genes were
# selected with AutoGeneS and saved, before running the Statescope framework and saving all steps.
# Key package versions:
#
# scRNAseq   v2.23.0

library(StatescopeR)
library(scRNAseq)

## Load scRNAseq
scRNAseq <- scRNAseq::SegerstolpePancreasData()

## remove duplicates gene names
scRNAseq <- scRNAseq[!duplicated(rownames(scRNAseq)), ]

## Preprocess scRNAseq
scRNAseq$donor <- scRNAseq$individual
scRNAseq$label <- scRNAseq$`cell type`

## Subset to 2 healthy and  type 2 diabetes samples
scRNAseq = scRNAseq[,scRNAseq$donor %in% c('H2', 'H3',
                                           'T2D1', 'T2D2')]

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

save(pseudobulk, file = 'inst/extdata/example_pseudobulk.RData')

## Create signature from scRNAseq for deconvolution
signature <- create_signature(scRNAseq, hvg_genes = TRUE, n_hvg_genes = 30L)

save(signature, file = 'inst/extdata/example_signature.RData')

## Select genes optimized for deconvolution (small number of genes for speed)
selected_genes <- select_genes(scRNAseq, 15L, n_hvg_genes = 30L)

save(selected_genes, file = 'inst/extdata/example_selected_genes.RData')

## Optionally create prior expectation
prior <- gather_true_fractions(scRNAseq) # Use True sc fractions for this
prior[rownames(prior) != "ductal cell", ] <- NA # Keep only ductal cell

## Tranpose it to nSample x nCelltype
prior <- t(prior)

## Perform Deconvolution with BLADE, refine gene expression estimates
Statescope <- BLADE_deconvolution(
    signature, pseudobulk, selected_genes,
    prior, 2L, Nrep = 2L
)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Deconvolved.RData')

Statescope <- Refinement(Statescope, signature, pseudobulk, 2L)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Refined.RData')

## Discover states
Statescope <- StateDiscovery(Statescope, k=2L, Ncores = 2L)

## Save to RData
save(Statescope, file = 'inst/extdata/example_Statescope_Discovered.RData')

