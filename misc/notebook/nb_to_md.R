## Convert an executed .ipynb into GitHub-flavoured Markdown with figure files.
## GitHub's notebook viewer runs its math pass over code cells, so R code with `$` is mangled there; the
## Markdown rendering (code fences protect `$`) is what the README links to.
## usage: Rscript nb_to_md.R notebook.ipynb out.md figure_dir
args <- commandArgs(TRUE); nb.file <- args[1]; md.file <- args[2]; fig.dir <- args[3]
nb <- jsonlite::fromJSON(nb.file, simplifyVector = FALSE)
unlink(fig.dir, recursive = TRUE); dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
rel <- file.path(basename(fig.dir), "")
out <- character(); k <- 0L
for (cell in nb$cells) {
  src <- paste(unlist(cell$source), collapse = "")
  if (cell$cell_type == "markdown") { out <- c(out, gsub("\\\\\\$", "$", src), ""); next }
  out <- c(out, "```r", src, "```", "")
  for (o in cell$outputs) {
    if (o$output_type == "stream") {
      txt <- paste(unlist(o$text), collapse = "")
      out <- c(out, "```", sub("\n$", "", txt), "```", "")
    } else if (o$output_type == "display_data" && !is.null(o$data$`image/png`)) {
      k <- k + 1L; f <- sprintf("fig-%02d.png", k)
      writeBin(jsonlite::base64_dec(o$data$`image/png`), file.path(fig.dir, f))
      out <- c(out, sprintf("![](%s%s)", rel, f), "")
    } else if (o$output_type == "error") {
      out <- c(out, "```", paste0("Error: ", o$evalue), "```", "")
    }
  }
}
writeLines(out, md.file)
cat("written", md.file, "with", k, "figures in", fig.dir, "\n")
