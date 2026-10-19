#!/usr/bin/env Rscript
# sync-doi-metadata.R -- make post/<slug>/_metadata.yml agree with the master dois.tsv.
#
# dois.tsv is the single source of truth.  The per-post _metadata.yml files are derived
# copies (Quarto reads the DOI from there to show it on the page).  This script
#   * adds   doi: "<doi>"  where a post's _metadata.yml has none (creating the file if needed);
#   * NEVER overwrites a doi that is already there: a different value is reported as a conflict;
#   * reports any _metadata.yml that carries a doi for a slug that is not in dois.tsv.
#
#   Rscript sync-doi-metadata.R          check only: exit 1 if anything is missing or wrong
#   Rscript sync-doi-metadata.R --write  add the missing lines; exit 1 only on conflicts/strays
# Run from the project root.  Handy in _quarto.yml:  pre-render: Rscript sync-doi-metadata.R

args <- commandArgs(trailingOnly = TRUE)
write_mode <- "--write" %in% args
master <- "dois.tsv"
header <- "slug\tdoi\tcreated_utc\ttitle"
die  <- function(...) { cat("ERROR:", ..., "\n", file = stderr()); quit(save = "no", status = 1) }
info <- function(...) cat(..., "\n", sep = "")
if (!requireNamespace("yaml", quietly = TRUE)) die("install the yaml package")
if (!file.exists("_quarto.yml")) die("run this from the project root")
if (!file.exists(master)) die(master, " not found")

norm_doi <- function(x) tolower(sub("^https?://(dx\\.)?doi\\.org/", "", trimws(as.character(x))))
ln <- readLines(master, warn = FALSE, encoding = "UTF-8")
if (ln[1] != header) die("unexpected header in ", master)
f <- strsplit(ln[-1], "\t", fixed = TRUE)
m <- setNames(vapply(f, `[`, "", 2), vapply(f, `[`, "", 1))      # slug -> doi
if (anyDuplicated(names(m))) die("duplicate slug in ", master)

missing <- character(); conflicts <- character(); strays <- character()
for (s in names(m)) {
  d <- file.path("post", s)
  if (!dir.exists(d)) next                                        # orphan: make-dois.R reports those
  p <- file.path(d, "_metadata.yml")
  cur <- NULL
  if (file.exists(p)) {
    y <- tryCatch(yaml::yaml.load_file(p), error = function(e) die("cannot parse ", p, ": ", conditionMessage(e)))
    if (!is.null(y$doi)) cur <- norm_doi(y$doi)
  }
  if (is.null(cur))            missing   <- c(missing, s)
  else if (cur != m[[s]])      conflicts <- c(conflicts, sprintf("%s: file has %s, master has %s", s, cur, m[[s]]))
}
# _metadata.yml files carrying a DOI for a slug that is not in the master
for (p in Sys.glob("post/*/_metadata.yml")) {
  s <- basename(dirname(p))
  y <- tryCatch(yaml::yaml.load_file(p), error = function(e) NULL)
  if (!is.null(y$doi) && !(s %in% names(m))) strays <- c(strays, sprintf("%s: %s", s, norm_doi(y$doi)))
}

if (length(conflicts)) { info("CONFLICTS (not touched; resolve by hand, master wins):"); info("  ", conflicts, sep = "\n  ") }
if (length(strays))    { info("STRAY DOIs not in ", master, ":"); info("  ", strays, sep = "\n  ") }

if (length(missing)) {
  if (!write_mode) {
    info(length(missing), " post(s) lack a doi in _metadata.yml: ", paste(head(missing, 8), collapse = ", "),
         if (length(missing) > 8) ", ..." else "")
  } else {
    for (s in missing) {
      p <- file.path("post", s, "_metadata.yml")
      line <- sprintf('doi: "%s"\n', m[[s]])
      if (file.exists(p)) {
        raw <- readBin(p, "raw", n = file.size(p))
        if (length(raw) && raw[length(raw)] != as.raw(0x0a)) line <- paste0("\n", line)
      }
      cat(line, file = p, append = TRUE)
    }
    info("Added doi to ", length(missing), " _metadata.yml file(s).")
  }
}
bad <- length(conflicts) + length(strays) > 0
pending <- !write_mode && length(missing) > 0
if (!bad && !pending) info("OK: ", sum(dir.exists(file.path("post", names(m)))), " post(s) consistent with ", master)
quit(save = "no", status = if (bad || pending) 1 else 0)
