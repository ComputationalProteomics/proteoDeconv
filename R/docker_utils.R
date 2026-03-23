docker_is_available <- function() {
  status <- suppressWarnings(
    system2("docker", "--version", stdout = FALSE, stderr = FALSE)
  )

  identical(as.integer(status), 0L)
}

docker_host_path <- function(path) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  gsub("\\\\", "/", normalized_path)
}

docker_container_path <- function(host_path, container_dir) {
  paste0(container_dir, "/", basename(host_path))
}

redact_docker_args <- function(args) {
  value_positions <- which(args %in% c("--username", "--token")) + 1L
  value_positions <- value_positions[value_positions <= length(args)]
  args[value_positions] <- "<redacted>"
  args
}

redact_docker_output <- function(text) {
  secrets <- Sys.getenv(c("CIBERSORTX_EMAIL", "CIBERSORTX_TOKEN"), unset = "")

  for (secret in secrets[nzchar(secrets)]) {
    text <- gsub(secret, "<redacted>", text, fixed = TRUE)
  }

  text
}

docker_mount_spec <- function(host_dir, container_dir) {
  selinux_suffix <- if (identical(Sys.info()[["sysname"]], "Linux")) ":z" else ""

  paste0(docker_host_path(host_dir), ":", container_dir, selinux_suffix)
}

build_docker_run_args <- function(image, mounts = character(), args = character()) {
  mount_args <- unlist(
    lapply(mounts, function(mount) c("-v", mount)),
    use.names = FALSE
  )

  c("run", "--rm", mount_args, image, args)
}

run_docker_command <- function(args, use_sudo = FALSE, verbose = FALSE) {
  command <- if (use_sudo) "sudo" else "docker"
  system_args <- if (use_sudo) c("docker", args) else args
  stdout_file <- tempfile()
  stderr_file <- tempfile()

  on.exit(unlink(c(stdout_file, stderr_file), force = TRUE), add = TRUE)

  if (verbose) {
    message(
      "Docker command:\n",
      paste(c(command, redact_docker_args(system_args)), collapse = " "),
      "\n"
    )
  }

  command_output <- system2(
    command,
    system_args,
    stdout = stdout_file,
    stderr = stderr_file
  )

  if (verbose || !identical(as.integer(command_output), 0L)) {
    output <- c(
      readLines(stdout_file, warn = FALSE),
      readLines(stderr_file, warn = FALSE)
    )

    if (length(output) > 0) {
      message(paste(redact_docker_output(output), collapse = "\n"))
    }
  }

  command_output
}
