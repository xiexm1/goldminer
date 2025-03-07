#!/usr/bin/env Rscript
args <- commandArgs(T)
prefix <- args[1]
out_dir <- args[2]

if (!require("tidyr")) {
    stop("Error: package 'tidyr' is not installed. Please install it first.")
}

clu <- read.table(paste0(out_dir, "/", prefix, ".csv"),
    col.names = c("spec", "cluid", "clu"), fill = T
)

clu$cluid <- paste0("clu:", clu$cluid)

clu_w <- tidyr::pivot_wider(clu,
    names_from = "spec",
    values_from = "cluid",
    values_fn = list(cluid = list),
    values_fill = list(cluid = list("."))
)

data.table::fwrite(clu_w, file = paste0(out_dir, "/", prefix, ".matrix"), sep = "	", row.names = F)