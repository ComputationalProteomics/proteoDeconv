gene_symbols_map <- HGNChelper::getCurrentHumanMap()

usethis::use_data(gene_symbols_map, overwrite = TRUE)
