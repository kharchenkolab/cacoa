#!/usr/bin/env bash
# Render the executed walkthrough as GitHub-flavoured Markdown: vignettes/walkthrough_short.md + walkthrough_short_files/.
# The Rmd is the single source; its chunks run only when the simulated dataset is present (see the setup chunk).
# usage: [R_LIBS=<lib>] misc/notebook/render_walkthrough.sh [path/to/shifts_sim_objects.rds] [n.cores]
# Needs pandoc >= 2.11.2 for rmarkdown::github_document; the system pandoc 2.5 is too old, so the RStudio/quarto
# bundled one is used when present.
set -euo pipefail
cd "$(dirname "$0")/../.."
sim=${1:-../test/shifts_sim_objects.rds}; cores=${2:-16}
[[ -f "$sim" ]] || { echo "simulated dataset not found: $sim" >&2; exit 1; }
sim=$(readlink -f "$sim")
for p in "${RSTUDIO_PANDOC:-}" /usr/lib/rstudio-server/bin/quarto/bin/tools /opt/quarto/bin/tools; do
  [[ -n "$p" && -x "$p/pandoc" ]] && { export RSTUDIO_PANDOC="$p"; break; }
done
Rscript -e "rmarkdown::render('vignettes/walkthrough_short.Rmd', output_format = rmarkdown::github_document(toc = TRUE, html_preview = FALSE),
  output_file = 'walkthrough_short.md', params = list(sim.file = '$sim', n.cores = ${cores}L), envir = new.env())"
