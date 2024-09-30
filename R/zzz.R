.onLoad <- function(libname, pkgname) {
  # if (Sys.getenv("PROTEODECONV_CACHE") == "") {
  #   db <- memoise::cache_filesystem(tools::R_user_dir("proteoDeconv", "cache"))
  # } else {
  #   db <- memoise::cache_filesystem(Sys.getenv("PROTEODECONV_CACHE"))
  # }
  # deconvolute <<- memoise::memoise(deconvolute, cache = db)
  # handle_missing_values <<- memoise::memoise(handle_missing_values, cache = db)
  # create_signature_matrix <<- memoise::memoise(create_signature_matrix, cache = db)
}
