suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out <- file.path(root, "Reference_Data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

terms <- c(
  "GO:0001501" = "skeletal_system_development",
  "GO:0001503" = "ossification",
  "GO:0045444" = "fat_cell_differentiation",
  "GO:0048762" = "mesenchymal_cell_differentiation"
)

go_map <- as.list(org.Hs.egGO2ALLEGS)
for (go_id in names(terms)) {
  entrez <- unique(go_map[[go_id]])
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = entrez,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  symbols <- sort(unique(na.omit(unname(symbols))))
  writeLines(symbols, file.path(out, paste0(go_id, "_", terms[[go_id]], ".txt")))
  cat(go_id, terms[[go_id]], length(symbols), "genes\n")
}
