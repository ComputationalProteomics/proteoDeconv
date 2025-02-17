#' Create Signature Matrix using CIBERSORTx
#'
#' This function generates a signature matrix from a reference sample and cell type classes using the CIBERSORTx Docker tool.
#'
#' @param refsample A data frame containing the reference profiles (with replicates). Required.
#' @param phenoclasses A data frame containing the cell type classes. Required.
#' @param g_min Minimum number of genes per cell type in the signature matrix. Default is 200.
#' @param g_max Maximum number of genes per cell type in the signature matrix. Default is 400.
#' @param q_value Q-value threshold for differential expression. Default is 0.01.
#' @param replicates Number of replicates to use for building the reference file. Default is 5.
#' @param sampling Fraction of available gene expression profiles selected by random sampling. Default is 0.5.
#' @param fraction Fraction of cells of the same identity showing evidence of expression. Default is 0.75.
#' @param filter Logical indicating whether to remove non-hematopoietic genes. Default is FALSE.
#' @param verbose Logical indicating whether to print detailed output. Default is FALSE.
#' @param QN Logical indicating whether to run quantile normalization. Default is FALSE.
#' @param single_cell Logical indicating whether to create signature from scRNA-Seq data. Default is FALSE.
#' @param ... Additional arguments passed to the function.
#'
#' @return A data frame containing the generated signature matrix.
#'
#' @details This function uses the CIBERSORTx Docker image to construct a signature matrix based on the input reference profiles and cell type classifications.
#'
#' @importFrom readr write_tsv read_tsv
#' @importFrom glue glue
#' @importFrom withr with_tempdir
#' @export
create_signature_matrix <- function(
    refsample = NULL, phenoclasses = NULL,
    g_min = 200,
    g_max = 400,
    q_value = 0.01,
    replicates = 5,
    sampling = 0.5,
    fraction = 0.75,
    filter = FALSE,
    verbose = FALSE,
    QN = FALSE,
    single_cell = FALSE,
    ...) {
  username <- Sys.getenv("CIBERSORTX_EMAIL")
  token <- Sys.getenv("CIBERSORTX_TOKEN")
  if (username == "" || token == "") {
    stop("CIBERSORTX_EMAIL and CIBERSORTX_TOKEN environment variables must be set.")
  }
  if (system("docker --version", intern = TRUE, ignore.stderr = TRUE) == 127) {
    stop("Docker is not installed or not running.")
  }

  withr::with_tempdir({
    input_dir <- tempdir()
    output_dir <- tempdir()

    refsample_file <- if (!is.null(refsample)) tempfile(tmpdir = input_dir, fileext = ".txt") else NULL
    phenoclasses_file <- if (!is.null(phenoclasses)) tempfile(tmpdir = input_dir, fileext = ".txt") else NULL

    if (!is.null(refsample)) readr::write_tsv(refsample, refsample_file)
    if (!is.null(phenoclasses)) readr::write_tsv(phenoclasses, phenoclasses_file, col_names = FALSE)

    docker_command <- glue::glue(
      "{if (Sys.info()[['sysname']] == 'Linux') 'sudo ' else ''}docker run --rm -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
      "--verbose {ifelse(verbose, 'TRUE', 'FALSE')} ",
      "--username {username} ",
      "--token {token} ",
      "{if (!is.null(refsample_file)) paste('--refsample', refsample_file) else ''} ",
      "{if (!is.null(phenoclasses_file)) paste('--phenoclasses', phenoclasses_file) else ''} ",
      "--G.min {g_min} ",
      "--G.max {g_max} ",
      "--q.value {q_value} ",
      "--filter {ifelse(filter, 'TRUE', 'FALSE')} ",
      "--QN {ifelse(QN, 'TRUE', 'FALSE')} ",
      "--k.max 999 ",
      "--remake FALSE ",
      "--replicates {replicates} ",
      "--sampling {sampling} ",
      "--fraction {fraction} ",
      "{ifelse(single_cell, '--single_cell TRUE', '')}"
    )
    if (verbose) {
      cat("Docker command:\n", docker_command, "\n")
    }
    command_output <- system(docker_command)
    if (command_output != 0) {
      stop(glue::glue("CIBERSORTx failed. Error code: {command_output}."))
    }
    refsample_id <- sub("\\.txt$", "", basename(refsample_file))
    signature_matrix_file <- grep(paste0("CIBERSORTx_", refsample_id),
      list.files(output_dir, pattern = "\\.txt$", full.names = TRUE),
      value = TRUE, fixed = TRUE
    )[1]
    if (is.na(signature_matrix_file)) {
      stop("No signature matrix file found in the output directory.")
    }
    signature_matrix <- readr::read_tsv(signature_matrix_file, show_col_types = FALSE)
  })
}

#' Create Phenoclasses for Immune Cells
#'
#' Generates a phenotype classification matrix from immune cell data,
#' for use with the CIBERSORTx signature generation function.
#' Each cell in the matrix is encoded as follows:
#' 0 for 'Unknown' mappings,
#' 1 if the cell type in the column matches the current group, or
#' 2 if it does not match.
#'
#' @param immune_cells A data frame of immune cell profiles. Must include a "Genes" column.
#' @param mapping_rules An optional named list of regex patterns for grouping; passed to map_cell_groups.
#' @param verbose Logical. If TRUE, displays additional messages during processing.
#'
#' @return A tibble with a "cell_type" column (the groups) and columns corresponding to the original profiles.
#'
#' @export
create_phenoclasses <- function(immune_cells,
                                mapping_rules = NULL,
                                gene_column = NULL,
                                verbose = FALSE) {
  immune_cells <- handle_input_data(immune_cells, gene_column = gene_column)

  data <- immune_cells %>% select(-Genes)
  cols <- colnames(data)

  group_mapped <- map_cell_groups(cols, mapping_rules, default_group = "Unknown", verbose = verbose)
  cell_type_to_group <- setNames(group_mapped, cols)


  valid_cols <- names(cell_type_to_group)
  valid_groups <- unique(cell_type_to_group)

  valid_groups <- setdiff(valid_groups, "Unknown")

  phenotype_classes <- matrix(NA,
    nrow = length(valid_groups),
    ncol = length(valid_cols),
    dimnames = list(valid_groups, valid_cols)
  )

  for (group in valid_groups) {
    for (j in seq_along(valid_cols)) {
      current_group <- cell_type_to_group[[valid_cols[j]]]
      if (current_group == "Unknown") {
        phenotype_classes[group, j] <- 0
      } else if (current_group == group) {
        phenotype_classes[group, j] <- 1
      } else {
        phenotype_classes[group, j] <- 2
      }
    }
  }

  phenotype_classes_df <- as.data.frame(phenotype_classes) %>%
    rownames_to_column(var = "cell_type") %>%
    as_tibble()

  return(phenotype_classes_df)
}
