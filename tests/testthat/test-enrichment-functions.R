
test_that("contingencyTable works", {

    se <- createSE()

    expect_s4_class(se, "SummarizedExperiment")
})

cols = "up"
for (i in seq_along(cols)) {
    row_data[cols[i]] <- factor(row_data[[cols[i]]], levels = c(TRUE, FALSE), labels = c(cols[i], paste0("not_", cols[i])))
}
