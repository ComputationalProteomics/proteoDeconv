test_that("docker helpers normalize host paths for bind mounts", {
  docker_host_path <- getFromNamespace("docker_host_path", "proteoDeconv")
  docker_mount_spec <- getFromNamespace("docker_mount_spec", "proteoDeconv")

  host_path <- if (.Platform$OS.type == "windows") {
    "C:\\temp\\proteoDeconv\\RtmpABC123"
  } else {
    "/tmp/Rtmp84NBvE"
  }

  normalized_host_path <- docker_host_path(host_path)
  expect_false(grepl("\\\\", normalized_host_path))

  mount_spec <- docker_mount_spec(host_path, "/src/data")
  expect_match(mount_spec, ":/src/data")

  if (identical(Sys.info()[["sysname"]], "Linux")) {
    expect_true(grepl(":z$", mount_spec))
  } else {
    expect_false(grepl(":z$", mount_spec))
  }
})

test_that("docker helpers redact sensitive arguments in logged commands", {
  redact_docker_args <- getFromNamespace("redact_docker_args", "proteoDeconv")

  args <- c(
    "run", "--username", "author@example.org", "--token", "super-secret-token"
  )

  expect_equal(
    redact_docker_args(args),
    c("run", "--username", "<redacted>", "--token", "<redacted>")
  )
})

test_that("create_signature_matrix passes container-visible file paths to Docker", {
  captured_args <- NULL
  tempfile_counter <- 0L

  withr::local_envvar(c(
    CIBERSORTX_EMAIL = "user@example.com",
    CIBERSORTX_TOKEN = "token"
  ))

  testthat::local_mocked_bindings(
    system2 = function(command, args, ...) {
      if (identical(command, "docker") && identical(args, "--version")) {
        return(0L)
      }

      captured_args <<- args
      0L
    },
    tempfile = function(pattern = "file", tmpdir = tempdir(), fileext = "") {
      tempfile_counter <<- tempfile_counter + 1L
      file.path(tmpdir, paste0("file", tempfile_counter, fileext))
    },
    list.files = function(path, pattern, full.names, ...) {
      file.path(path, "CIBERSORTx_file2.txt")
    },
    .package = "base"
  )

  testthat::local_mocked_bindings(
    read_tsv = function(file, show_col_types = FALSE, ...) {
      tibble::tibble(GeneSymbol = "Gene1", CellType = 1)
    },
    .package = "readr"
  )

  refsample <- matrix(
    1,
    nrow = 1,
    ncol = 1,
    dimnames = list("Gene1", "Sample1")
  )
  phenoclasses <- matrix(
    1,
    nrow = 1,
    ncol = 1,
    dimnames = list("CellType", "Sample1")
  )

  result <- create_signature_matrix(
    refsample = refsample,
    phenoclasses = phenoclasses
  )

  expect_true(is.matrix(result))

  refsample_arg <- captured_args[match("--refsample", captured_args) + 1]
  phenoclasses_arg <- captured_args[match("--phenoclasses", captured_args) + 1]

  expect_match(refsample_arg, "^/src/data/")
  expect_match(phenoclasses_arg, "^/src/data/")
  expect_false(grepl(tempdir(), refsample_arg, fixed = TRUE))
  expect_false(grepl(tempdir(), phenoclasses_arg, fixed = TRUE))
})

test_that("deconvolute_cibersortx passes container-visible file paths to Docker", {
  captured_args <- NULL
  tempfile_counter <- 0L

  withr::local_envvar(c(
    CIBERSORTX_EMAIL = "user@example.com",
    CIBERSORTX_TOKEN = "token"
  ))

  testthat::local_mocked_bindings(
    system2 = function(command, args, ...) {
      if (identical(command, "docker") && identical(args, "--version")) {
        return(0L)
      }

      captured_args <<- args
      0L
    },
    tempfile = function(pattern = "file", tmpdir = tempdir(), fileext = "") {
      tempfile_counter <<- tempfile_counter + 1L
      file.path(tmpdir, paste0("file", tempfile_counter, fileext))
    },
    file.exists = function(path) {
      grepl("CIBERSORTx_.*_Results\\.txt$", path)
    },
    .package = "base"
  )

  testthat::local_mocked_bindings(
    read_tsv = function(file, show_col_types = FALSE, ...) {
      tibble::tibble(Mixture = "Sample1", CellType = 0.75)
    },
    .package = "readr"
  )

  data <- matrix(
    1,
    nrow = 1,
    ncol = 1,
    dimnames = list("Gene1", "Sample1")
  )
  signature <- matrix(
    1,
    nrow = 1,
    ncol = 1,
    dimnames = list("Gene1", "CellType")
  )

  result <- deconvolute_cibersortx(data = data, signature = signature)

  expect_true(is.matrix(result))

  mixture_arg <- captured_args[match("--mixture", captured_args) + 1]
  signature_arg <- captured_args[match("--sigmatrix", captured_args) + 1]
  source_geps_arg <- captured_args[match("--sourceGEPs", captured_args) + 1]

  expect_match(mixture_arg, "^/src/data/")
  expect_match(signature_arg, "^/src/data/")
  expect_match(source_geps_arg, "^/src/data/")
  expect_false(grepl(tempdir(), mixture_arg, fixed = TRUE))
  expect_false(grepl(tempdir(), signature_arg, fixed = TRUE))
})
