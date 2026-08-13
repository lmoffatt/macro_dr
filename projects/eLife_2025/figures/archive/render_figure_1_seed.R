#!/usr/bin/env Rscript
# Render Figure 1 (and its full-recording companion) from one seed of the sweep produced by
# ops/local/sweep_figure_1_seeds.sh, without touching figure_1.Rmd or figure_1_panels.R.
#
# figure_1_panels.R reads its dumps as "../data/figure_1_likelihood_diagnostic_<ALGO>.csv",
# relative to the working directory. So rendering a given seed is a matter of running from a
# directory whose parent holds that seed's data/: each seed dir already has figures/data/, and
# this script creates figures/render/ beside it and works there. Nothing permanent changes,
# and the committed figure_1.Rmd keeps rendering the canonical data as before.
#
# Usage:
#   Rscript render_figure_1_seed.R 7
#   Rscript render_figure_1_seed.R 7 --sel 1,2
#   Rscript render_figure_1_seed.R 7 --cols LSE,NR,R,IR --install
#
#   --sel      comma-separated sample indices for the zoom panel (default 3,4)
#   --full     comma-separated sample indices for the companion (default 0,1,2,3,4)
#   --cols     roster, left to right (default LSE,NR,R,IR)
#   --outroot  sweep root (default ../../ops/local/figure_1_seeds)
#   --install  also copy the two PDFs over Figure_1.pdf / Figure_S1.pdf in this directory

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: render_figure_1_seed.R <seed> [--sel 3,4] [--cols ...] [--install]")

seed <- as.integer(args[1])
if (is.na(seed)) stop("first argument must be the seed number")

opt <- function(name, default) {
  k <- match(paste0("--", name), args)
  if (is.na(k) || k == length(args)) default else args[k + 1]
}
has <- function(name) paste0("--", name) %in% args
as_idx <- function(s) as.integer(strsplit(s, ",", fixed = TRUE)[[1]])

sel     <- as_idx(opt("sel", "3,4"))
full    <- as_idx(opt("full", "0,1,2,3,4"))
cols    <- strsplit(opt("cols", "LSE,NR,R,IR"), ",", fixed = TRUE)[[1]]
outroot <- opt("outroot", "../../ops/local/figure_1_seeds")

# Absolute paths first: everything below runs from a different working directory.
here    <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
panels  <- file.path(here, "figure_1_panels.R")
outroot <- normalizePath(file.path(here, outroot), mustWork = FALSE)
if (!dir.exists(outroot)) outroot <- normalizePath(opt("outroot", outroot), mustWork = TRUE)
seeddir <- file.path(outroot, paste0("seed_", seed))
datadir <- file.path(seeddir, "figures", "data")
stopifnot(file.exists(panels), dir.exists(datadir))

rdir <- file.path(seeddir, "figures", "render")
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)

# Prefixed, because source()ing the panels script brings its own globals into this
# environment and some of them are short names (fs, d, noise, cols would all collide).
.rf1_zoom <- paste0("Figure_1_seed", seed, ".pdf")
.rf1_full <- paste0("Figure_S1_seed", seed, ".pdf")
.rf1_sel  <- sel
.rf1_win  <- full
.rf1_cols <- cols
.rf1_rdir <- rdir
.rf1_here <- here
.rf1_seed <- seed
.rf1_install <- has("install")

.rf1_owd <- setwd(.rf1_rdir)
on.exit(setwd(.rf1_owd), add = TRUE)

message("[render] seed ", .rf1_seed, "  roster ",
        paste(.rf1_cols, collapse = "/"),
        "  zoom ", paste(.rf1_sel, collapse = ","))

# Brings its own globals into this environment; everything used after this point
# is prefixed so that fs, d, noise and friends cannot shadow it.
source(panels, local = FALSE)

build_figure(.rf1_sel, .rf1_zoom, .rf1_cols)
build_figure(.rf1_win, .rf1_full, .rf1_cols)

setwd(.rf1_owd)

cat("\n[render] written:\n")
for (.rf1_p in file.path(.rf1_rdir, c(.rf1_zoom, .rf1_full))) {
  if (file.exists(.rf1_p)) cat(sprintf("  %-28s %s\n", basename(.rf1_p), .rf1_p))
}

if (.rf1_install) {
  file.copy(file.path(.rf1_rdir, .rf1_zoom),
            file.path(.rf1_here, "Figure_1.pdf"), overwrite = TRUE)
  file.copy(file.path(.rf1_rdir, .rf1_full),
            file.path(.rf1_here, "Figure_S1.pdf"), overwrite = TRUE)
  cat("[render] installed as Figure_1.pdf and Figure_S1.pdf (seed ", .rf1_seed,
      " — record it in the caption)\n", sep = "")
}
