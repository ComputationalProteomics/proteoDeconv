# run inside the manuscript targets pipeline
# todo: make this independent of the pipeline

# library(targets)
# library(tidyverse)
# tar_load_globals()

# mix <- tar_read(dia_pg) |>
#     read_tsv() |>
#     select(contains(c("Genes", "Mix")))

# cd8_mono <- tar_read(dia_pg) |>
#     read_tsv() |>
#     select(contains(c("Genes", "CD8", "Mono")))

# saveRDS(mixes_data, "mixed_samples_matrix.rds")
# saveRDS(cd8_mono, "cd8_mono.rds")
