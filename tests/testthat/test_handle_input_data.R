test_that("SummarizedExperiment input is handled correctly", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(data = matrix(1:4, nrow = 2, dimnames = list(c("GeneA", "GeneB"), NULL)))
    )
    result <- handle_input_data(se)
    expect_true("Genes" %in% colnames(result))
    expect_equal(nrow(result), 2)
})

test_that("data.frame input is handled correctly", {
    df <- data.frame(A = 1:2, B = 3:4, row.names = c("GeneA", "GeneB"))
    result <- handle_input_data(df)
    expect_true("Genes" %in% colnames(result))
    expect_equal(nrow(result), 2)
})

test_that("invalid input throws an error", {
    expect_error(handle_input_data(list(a = 1, b = 2)))
})
