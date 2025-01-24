test_that("update_gene_symbols handles data frame input", {
    df <- data.frame(Genes = c("GeneA", "GeneB"))
    result <- update_gene_symbols(df)
    expect_s3_class(result, "data.frame")
})

test_that("update_gene_symbols handles vector input", {
    genes <- c("GeneA", "GeneB")
    result <- update_gene_symbols(genes)
    expect_vector(result)
})

# test_that("update_gene_symbols fails on invalid input", {
#     expect_error(update_gene_symbols(list("not", "supported")))
# })
