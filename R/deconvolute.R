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
#'   \item Other algorithms supported by the `immunedeconv` package.
#' }
#'
#' @importFrom tibble column_to_rownames
#' @importFrom immunedeconv deconvolute
#' @export
deconvolute <- function(algorithm, data, signature_matrix = NULL, ...) {
    if (algorithm == "cibersortx") {
        return(deconvolute_cibersortx(data, signature_matrix, ...))
    }

    if (algorithm == "cibersort") {
        return(deconvolute_cibersort(data, signature_matrix))
    }

    data <- data |> tibble::column_to_rownames("Genes")

    if (algorithm == "epic") {
        return(deconvolute_epic(data, signature_matrix, ...))
    }

    result <- immunedeconv::deconvolute(data, algorithm)
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
#'
#' @return A data frame containing the deconvolution results.
#'
#' @details
#' This function requires the CIBERSORTx Docker image to be installed. It also requires the user to set the \code{CIBERSORTX_EMAIL} and \code{CIBERSORTX_TOKEN} environment variables with their CIBERSORTx credentials.
#'
#' @export
deconvolute_cibersortx <- function(data, signature_matrix, perm = 1, rmbatch_S_mode = FALSE, source_GEPs = NULL, use_cibersortx = TRUE, rmbatch_B_mode = FALSE, QN = FALSE, absolute = FALSE, abs_method = "sig.score") {
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
            "docker run -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
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

deconvolute_cibersort <- function(data, signature_matrix, ...) {
    cibersort_result <- immunedeconv::deconvolute_cibersort_custom(
        data |> tibble::column_to_rownames(colnames(data)[1]),
        signature_matrix |> tibble::column_to_rownames(colnames(signature_matrix)[1])
    )
    tibble::as_tibble(cibersort_result |> t(), rownames = "Mixture") |> convert_cibersortx_output()
}

deconvolute_epic <- function(data, signature_matrix, signature_matrix_variance, signature_genes, default_reference = EPIC::TRef, ...) {
    result <- immunedeconv::deconvolute_epic_custom(
        data,
        signature_matrix = signature_matrix,
        genes_var = signature_matrix_variance,
        signature_genes = signature_genes
    )
    tibble::as_tibble(result, rownames = "cell_type")
}
