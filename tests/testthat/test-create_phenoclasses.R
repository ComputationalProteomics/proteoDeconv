test_that("create_phenoclasses correctly creates phenotype classification matrix", {
  gene_names <- c("Gene1", "Gene2", "Gene3")
  sample_names <- c("T_cell_CD4", "T_cell_CD8", "B_cell", "Macrophage_M1")

  immune_cells_mat <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    nrow = 3,
    ncol = 4,
    dimnames = list(gene_names, sample_names)
  )

  mapping_rules <- list(
    "T_cell" = c("T_cell", "CD4", "CD8"),
    "B_cell" = c("B_cell"),
    "Macrophage" = c("Macrophage", "M1")
  )

  result_mat <- create_phenoclasses(
    immune_cells_mat,
    mapping_rules = mapping_rules,
    return_format = "matrix"
  )

  expect_true(is.matrix(result_mat))
  expect_equal(nrow(result_mat), 3)
  expect_equal(ncol(result_mat), 4)
  expect_setequal(rownames(result_mat), c("T_cell", "B_cell", "Macrophage"))
  expect_setequal(
    colnames(result_mat),
    c("T_cell_CD4", "T_cell_CD8", "B_cell", "Macrophage_M1")
  )

  expect_equal(result_mat["T_cell", "T_cell_CD4"], 1)
  expect_equal(result_mat["T_cell", "T_cell_CD8"], 1)
  expect_equal(result_mat["T_cell", "B_cell"], 2)
  expect_equal(result_mat["T_cell", "Macrophage_M1"], 2)

  expect_equal(result_mat["B_cell", "T_cell_CD4"], 2)
  expect_equal(result_mat["B_cell", "T_cell_CD8"], 2)
  expect_equal(result_mat["B_cell", "B_cell"], 1)
  expect_equal(result_mat["B_cell", "Macrophage_M1"], 2)

  expect_equal(result_mat["Macrophage", "T_cell_CD4"], 2)
  expect_equal(result_mat["Macrophage", "T_cell_CD8"], 2)
  expect_equal(result_mat["Macrophage", "B_cell"], 2)
  expect_equal(result_mat["Macrophage", "Macrophage_M1"], 1)

  result_tibble <- create_phenoclasses(
    immune_cells_mat,
    mapping_rules = mapping_rules,
    return_format = "tibble"
  )

  expect_true(tibble::is_tibble(result_tibble))
  expect_equal(nrow(result_tibble), 3)
  expect_true("cell_type" %in% colnames(result_tibble))
  expect_setequal(result_tibble$cell_type, c("T_cell", "B_cell", "Macrophage"))
})


test_that("create_phenoclasses handles verbose output", {
  gene_names <- c("Gene1", "Gene2")
  sample_names <- c("T_cell", "Unknown_cell")

  immune_cells_mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    ncol = 2,
    dimnames = list(gene_names, sample_names)
  )

  mapping_rules <- list(
    "T_cell" = c("T_cell")
  )

  expect_no_error(
    create_phenoclasses(
      immune_cells_mat,
      mapping_rules = mapping_rules,
      verbose = TRUE
    )
  )

  expect_no_error(
    create_phenoclasses(
      immune_cells_mat,
      mapping_rules = mapping_rules,
      verbose = TRUE
    )
  )
})
