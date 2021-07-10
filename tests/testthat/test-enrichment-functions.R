
test_that("contingencyTable works", {

    se <- exampleSE()

    expect_s4_class(se, "SummarizedExperiment")
})
