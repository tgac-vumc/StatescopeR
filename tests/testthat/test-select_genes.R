## Check that input is a SingleCellExperiment
test_that("scRNAseq input is not a SingleCellExperiment", {
    expect_s4_class(data, 'SingleCellExperiment')
})

## Check that logcounts are actually log cp10k counts
test_that("logcounts are actually logged cp10k", {
    expect_equal(mean(colSums(exp(logcounts(data))-1)), 10000)
})
