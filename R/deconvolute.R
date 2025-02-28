#' Deconvolute Protein Data
#'
#' This function deconvolutes protein or gene expression data using various algorithms.
#'
#' @param algorithm A string specifying the deconvolution algorithm to use. Options are "cibersortx", "cibersort", "epic", and other algorithms supported by the `immunedeconv` package.
#' @param data A data frame containing protein or gene data. The data frame should have proteins/genes as rows and samples as columns.
#' @param signature_matrix An optional signature matrix required by some algorithms.
#' @param ... Additional arguments passed to the specific deconvolution functions.
#'
#' @return A data frame with the deconvolution results.
#'
#' @details The function supports multiple deconvolution algorithms:
#' \itemize{
#'   \item \strong{"cibersortx"}: Uses the `deconvolute_cibersortx` function.
#'   \item \strong{"cibersort"}: Uses the `deconvolute_cibersort` function.
#'   \item \strong{"epic"}: Uses the `deconvolute_epic` function.
#' }
#'
#' @importFrom tibble column_to_rownames
#' @export
deconvolute <- function(algorithm, data, signature_matrix = NULL, ...) {
    data <- handle_input_data(data, "Genes")

    if (algorithm == "cibersortx") {
        return(deconvolute_cibersortx(data, signature_matrix, ...))
    }

    if (algorithm == "cibersort") {
        return(deconvolute_cibersort(data, signature_matrix))
    }

    data <- data |> handle_input_data(gene_column = "Genes", as_tibble = FALSE)
}



#' Deconvolute using CIBERSORTx
#'
#' This function performs deconvolution using the CIBERSORTx Docker image.
#'
#' @param data A data frame containing the mixture data to be deconvoluted.
#' @param signature_matrix A data frame containing the signature matrix.
#' @param perm An integer specifying the number of permutations to be performed. Default is 1.
#' @param rmbatch_S_mode A logical value indicating whether to remove batch effects in source GEPs mode. Default is FALSE.
#' @param source_GEPs A data frame containing the source gene expression profiles. Required if \code{rmbatch_S_mode} is TRUE.
#' @param use_cibersortx A logical value indicating whether to use CIBERSORTx. Default is TRUE.
#' @param rmbatch_B_mode A logical value indicating whether to remove batch effects in bulk mode. Default is FALSE.
#' @param QN A logical value indicating whether to perform quantile normalization. Default is FALSE.
#' @param absolute A logical value indicating whether to use absolute mode. Default is FALSE.
#' @param abs_method A character string specifying the method to use for absolute mode. Default is "sig.score".
#' @param use_sudo A logical value indicating whether to use sudo for Docker commands. Default is FALSE.
#'
#' @return A data frame containing the deconvolution results.
#'
#' @details
#' This function requires the CIBERSORTx Docker image to be installed. It also requires the user to set the \code{CIBERSORTX_EMAIL} and \code{CIBERSORTX_TOKEN} environment variables with their CIBERSORTx credentials.
#'
#' @export
deconvolute_cibersortx <- function(data, signature_matrix, perm = 1, rmbatch_S_mode = FALSE, source_GEPs = NULL, use_cibersortx = TRUE, rmbatch_B_mode = FALSE, QN = FALSE, absolute = FALSE, abs_method = "sig.score", use_sudo = FALSE) {
    if (!use_cibersortx) {
        return(deconvolute_cibersort(data, signature_matrix))
    }

    if (absolute == TRUE) {
        stop("Absolute mode is not supported in the CIBERSORTx Docker image. Use the online version instead.")
    }

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

        input_data_file <- tempfile(tmpdir = input_dir)
        signature_matrix_file <- tempfile(tmpdir = input_dir)

        source_GEPs_file <- NULL
        if (rmbatch_S_mode) {
            source_GEPs_file <- tempfile(tmpdir = input_dir)
            readr::write_tsv(source_GEPs, source_GEPs_file)
        }

        readr::write_tsv(data, input_data_file)
        readr::write_tsv(signature_matrix, signature_matrix_file)

        label <- uuid::UUIDgenerate(TRUE)

        docker_command <- glue::glue(
            "{if (use_sudo) 'sudo ' else ''}docker run --rm -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
            "--verbose TRUE ",
            "--username {username} ",
            "--token {token} ",
            "--mixture {input_data_file} ",
            "--sigmatrix {signature_matrix_file} ",
            "--perm {perm} ",
            "--label {label} ",
            "--rmbatchBmode {ifelse(rmbatch_B_mode, 'TRUE', 'FALSE')} ",
            "--rmbatchSmode {ifelse(rmbatch_S_mode, 'TRUE', 'FALSE')} ",
            "--sourceGEPs {source_GEPs_file %||% signature_matrix_file} ",
            "--QN {ifelse(QN, 'TRUE', 'FALSE')} ",
            "--absolute {ifelse(absolute, 'TRUE', 'FALSE')} ",
            "--abs_method {abs_method} "
        )
        message("Docker command:\n", docker_command, "\n")

        command_output <- system(docker_command)

        if (command_output != 0) {
            stop(glue::glue("CIBERSORTx failed. Error code: {command_output}"))
        }

        result_file <- glue::glue("{output_dir}/CIBERSORTx_{label}_Results.txt")
        if (!file.exists(result_file)) {
            stop("No result file found in the output directory.")
        }

        cibersortx_result <- readr::read_tsv(result_file, show_col_types = FALSE)

        convert_cibersortx_output(cibersortx_result)
    })
}

#' Deconvolute Mixture Data using CIBERSORT
#'
#' Applies a custom CIBERSORT deconvolution on the input data using the provided signature matrix,
#' creating temporary input files and running the external CIBERSORT script in an isolated directory.
#' Results are returned as a tibble after removing unneeded columns and transposing.
#'
#' @param data A data frame containing mixture expression data.
#' @param signature_matrix A data frame with the signature matrix.
#' @param QN Logical indicating whether quantile normalization is performed (default FALSE).
#' @param absolute Logical indicating whether an absolute score is computed (default FALSE).
#' @param abs_method Method for absolute scoring if absolute is TRUE (default "sig.score").
#' @param ... Additional arguments passed to the CIBERSORT function.
#'
#' @return A tibble with deconvolution results.
#' @export
deconvolute_cibersort <- function(
    data,
    signature_matrix,
    QN = FALSE,
    absolute = FALSE,
    abs_method = "sig.score",
    ...) {
    data <- handle_input_data(data, as_tibble = FALSE)
    signature_matrix <- handle_input_data(signature_matrix, as_tibble = FALSE)

    if (!exists("CIBERSORT", mode = "function")) {
        stop("Function 'CIBERSORT' not found. Please load the CIBERSORT.R script.")
    }

    cibersort_result <- withr::with_tempdir({
        expr_file <- tempfile()
        sig_file <- tempfile()
        expr_tbl <- dplyr::as_tibble(data, rownames = "gene_symbol")
        sig_tbl <- dplyr::as_tibble(signature_matrix, rownames = "gene_symbol")
        readr::write_tsv(expr_tbl, file = expr_file)
        readr::write_tsv(sig_tbl, file = sig_file)

        extras <- rlang::dots_list(
            sig_file,
            expr_file,
            perm = 0,
            QN = QN,
            absolute = absolute,
            abs_method = abs_method,
            ...,
            .homonyms = "last"
        )
        cibersort_call <- rlang::call2(CIBERSORT, !!!extras)
        output <- eval(cibersort_call)

        pruned_output <- output %>%
            .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]
    })

    tibble::as_tibble(cibersort_result, rownames = "Mixture") |>
        convert_cibersortx_output()
}

convert_cibersortx_output <- function(data) {
    extra_cols <- c("P.value", "Correlation", "RMSE", "P-value")

    if (!("Mixture" %in% colnames(data))) {
        stop("No Mixture column found")
    }

    data |>
        dplyr::select(-tidyselect::any_of(extra_cols)) |>
        tidyr::pivot_longer(
            cols = -dplyr::all_of("Mixture"),
            names_to = "cell_type",
            values_to = "value"
        ) |>
        tidyr::pivot_wider(
            names_from = "Mixture",
            values_from = "value"
        )
}
