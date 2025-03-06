test_that("map_cell_groups correctly maps column names to cell groups", {
  column_names <- c(
    "T_cell_CD4",
    "B_cell",
    "Macrophage_M1",
    "Unknown_cell",
    "NK_cell"
  )

  mapping_rules <- list(
    "T_cell" = c("T_cell", "CD4", "CD8"),
    "B_cell" = c("B_cell", "B-cell"),
    "Macrophage" = c("Macrophage", "M1", "M2"),
    "NK_cell" = c("NK", "Natural Killer")
  )

  result <- map_cell_groups(
    column_names,
    mapping_rules,
    default_group = "Unknown"
  )

  expect_equal(result[1], "T_cell")
  expect_equal(result[2], "B_cell")
  expect_equal(result[3], "Macrophage")
  expect_equal(result[5], "NK_cell")

  expect_true(is.character(result[4]))
})

test_that("map_cell_groups handles case insensitivity", {
  column_names <- c("t_cell_cd4", "B_CELL", "Macrophage_m1")

  mapping_rules <- list(
    "T_cell" = c("T_cell", "CD4"),
    "B_cell" = c("B_cell"),
    "Macrophage" = c("Macrophage", "M1")
  )

  result <- map_cell_groups(column_names, mapping_rules)

  expected <- c("T_cell", "B_cell", "Macrophage")

  expect_equal(result, expected)
})

test_that("map_cell_groups returns default group for non-matches", {
  column_names <- c("Fibroblast", "Epithelial", "Endothelial")

  mapping_rules <- list(
    "T_cell" = c("T_cell", "CD4", "CD8"),
    "B_cell" = c("B_cell", "B-cell")
  )

  result <- map_cell_groups(
    column_names,
    mapping_rules,
    default_group = "Other"
  )

  expected <- c("Other", "Other", "Other")

  expect_equal(result, expected)
})

test_that("map_cell_groups throws error for non-character input", {
  expect_error(
    map_cell_groups(1:5, list("Group" = "pattern")),
    "`column_names` must be a character vector."
  )
})

test_that("map_cell_groups handles verbose output", {
  column_names <- c("T_cell_CD4", "Unknown_cell")

  mapping_rules <- list(
    "T_cell" = c("T_cell", "CD4")
  )

  expect_no_error(
    map_cell_groups(column_names, mapping_rules, verbose = TRUE)
  )
})
