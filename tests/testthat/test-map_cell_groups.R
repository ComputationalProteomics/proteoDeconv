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


test_that("map_cell_groups with use_regex=FALSE performs exact string matching", {
  column_names <- c(
    "T_cell_CD4",
    "B_cell",
    "Macrophage_M1",
    "NK_cell",
    "CD4"
  )

  mapping_rules <- list(
    "T_cell" = c("T_cell_CD4"),
    "B_cell" = c("B_cell"),
    "Macrophage" = c("Macrophage_M1"),
    "CD4_T_cell" = c("CD4")
  )

  result <- map_cell_groups(
    column_names,
    mapping_rules,
    use_regex = FALSE
  )

  expected <- c("T_cell", "B_cell", "Macrophage", "Unknown", "CD4_T_cell")

  expect_equal(result, expected)
})

test_that("map_cell_groups with use_regex=TRUE handles regex patterns", {
  column_names <- c(
    "T_cell_CD4",
    "B_cell",
    "Macrophage_M1",
    "NK_cell",
    "CD4"
  )

  mapping_rules <- list(
    "T_cell" = c("T_cell.*", "^CD4$"),
    "B_cell" = c("B[_-]cell"),
    "Macrophage" = c("Macrophage.*"),
    "NK_cell" = c("NK.*")
  )

  result <- map_cell_groups(
    column_names,
    mapping_rules,
    use_regex = TRUE
  )

  expected <- c("T_cell", "B_cell", "Macrophage", "NK_cell", "T_cell")

  expect_equal(result, expected)
})

test_that("map_cell_groups with use_regex=FALSE is case-sensitive", {
  column_names <- c("T_cell", "t_cell", "B_CELL", "b_cell")

  mapping_rules <- list(
    "T_cell" = c("T_cell"),
    "B_cell" = c("B_CELL")
  )

  result <- map_cell_groups(
    column_names,
    mapping_rules,
    use_regex = FALSE
  )

  expected <- c("T_cell", "Unknown", "B_cell", "Unknown")

  expect_equal(result, expected)
})

test_that("map_cell_groups with mixing use_regex modes gives different results", {
  column_names <- c("T_cell_CD4", "B-cell", "Mono")

  mapping_rules <- list(
    "T_cell" = c("T_cell.*"),
    "B_cell" = c("B-cell"),
    "Monocyte" = c("Mono.*")
  )

  result_regex <- map_cell_groups(
    column_names,
    mapping_rules,
    use_regex = TRUE
  )

  result_no_regex <- map_cell_groups(
    column_names,
    mapping_rules,
    use_regex = FALSE
  )

  expected_regex <- c("T_cell", "B_cell", "Monocyte")
  expected_no_regex <- c("Unknown", "B_cell", "Unknown")

  expect_equal(result_regex, expected_regex)
  expect_equal(result_no_regex, expected_no_regex)
  expect_false(identical(result_regex, result_no_regex))
})
