# test_that("duplication mode=slice works", {
#   data <- tibble::tribble(
#     ~Genes, ~sample1, ~sample2, ~sample3,
#     "A", 1, 2, 3,
#     "A", 2, 3, 3,
#     "C", 1, 2, 3
#   )
#   data_expected <- tibble::tribble(
#     ~Genes, ~sample1, ~sample2, ~sample3,
#     "A", 2, 3, 3,
#     "C", 1, 2, 3
#   )
#   expect_equal(
#     handle_duplicate_genes(data, duplicate_mode = "slice"),
#     data_expected
#   )
# })
