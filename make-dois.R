#!/usr/bin/env Rscript
# make-dois.R -- give each blog post a permanent Rogue Scholar DOI, exactly once.
#
# Master file: dois.tsv in the project root, columns: slug, doi, created_utc, title.
# It is APPEND-ONLY.  This script only ever adds lines at the end of the file; it
# never edits, reorders, or deletes an existing line, and it never generates a
# second DOI for a slug that already has one.
#
#   Rscript make-dois.R                  dry run (default): report only, create nothing
#   Rscript make-dois.R --write          append DOIs for posts that have none yet
#   Rscript make-dois.R --allow-orphans  proceed even when dois.tsv lists slugs that
#                                        no longer have a post folder (see below)
#   Rscript make-dois.R --prefix=10.59350   DOI prefix for NEW DOIs (this is the default)
#
# Environment: COMMONMETA=/path/to/commonmeta   (default: "commonmeta" on the PATH)
# Run from the project root (the folder containing _quarto.yml).

args <- commandArgs(trailingOnly = TRUE)
write_mode    <- "--write" %in% args
allow_orphans <- "--allow-orphans" %in% args
prefix <- sub("^--prefix=", "", grep("^--prefix=", args, value = TRUE))
if (length(prefix) != 1) prefix <- "10.59350"
cm     <- Sys.getenv("COMMONMETA", "commonmeta")
master <- "dois.tsv"
header <- "slug\tdoi\tcreated_utc\ttitle"

die  <- function(...) { cat("ERROR:", ..., "\n", file = stderr()); quit(save = "no", status = 1) }
warn <- function(...) cat("WARNING:", ..., "\n", file = stderr())
info <- function(...) cat(..., "\n", sep = "")

if (!requireNamespace("yaml", quietly = TRUE)) die("install the yaml package: install.packages('yaml')")
if (!file.exists("_quarto.yml")) die("run this from the project root (no _quarto.yml in ", getwd(), ")")

doi_any <- "^10\\.[0-9]{4,9}/[0-9a-z]{5}-[0-9a-z]{3}[0-9]{2}$"       # any prefix: shape of a Rogue Scholar DOI
doi_new <- paste0("^", gsub(".", "\\.", prefix, fixed = TRUE), "/[0-9a-z]{5}-[0-9a-z]{3}[0-9]{2}$")
slug_ok <- "^[A-Za-z0-9._-]+$"

norm_doi <- function(x) tolower(sub("^https?://(dx\\.)?doi\\.org/", "", trimws(x)))

# ---- commonmeta -------------------------------------------------------------
run_cm <- function(a) {
  tryCatch(suppressWarnings(system2(cm, a, stdout = TRUE, stderr = TRUE)),
           error = function(e) structure(character(), status = 127L))
}
# `commonmeta decode` prints the integer behind a valid DOI.  On a bad checksum the
# Python version prints 0 and still exits 0, so "valid" means: exit 0 AND a positive integer.
decode_ok <- function(doi) {
  out <- run_cm(c("decode", doi)); st <- attr(out, "status")
  if (!is.null(st) && st != 0) return(FALSE)
  out <- trimws(out); out <- out[nzchar(out)]
  length(out) == 1 && grepl("^[0-9]+$", out) && as.numeric(out) > 0
}
cm_selftest <- function() {
  # Example from the Rogue Scholar documentation (decodes to 227095303906); last digit changed = invalid.
  decode_ok("10.53731/6kfyy-nq280") && !decode_ok("10.53731/6kfyy-nq281")
}
cm_present <- nzchar(Sys.which(cm)) || file.exists(cm)
if (cm_present && !cm_selftest()) die("commonmeta self-test failed; refusing to continue (is '", cm, "' the right tool?)")
if (!cm_present) {
  if (write_mode) die("commonmeta not found ('", cm, "'). pip install commonmeta-py, or set COMMONMETA=/path")
  warn("commonmeta not found; a --write run would fail.")
}

# ---- read and validate the master file --------------------------------------
read_master <- function() {
  empty <- data.frame(slug = character(), doi = character(), created_utc = character(),
                      title = character(), stringsAsFactors = FALSE)
  if (!file.exists(master)) return(empty)
  raw <- readBin(master, "raw", n = file.size(master))
  if (!length(raw)) die(master, " exists but is empty; delete it only if you have never generated a DOI")
  if (raw[length(raw)] != as.raw(0x0a)) die(master, " does not end with a newline; fix by hand before appending")
  if (any(raw == as.raw(0x0d))) die(master, " contains carriage returns (CRLF); fix by hand")
  ln <- strsplit(rawToChar(raw), "\n", fixed = TRUE)[[1]]
  if (ln[1] != header) die("first line of ", master, " is not the expected header")
  ln <- ln[-1]
  if (any(!nzchar(ln))) die(master, " contains blank lines")
  f <- strsplit(ln, "\t", fixed = TRUE)
  if (any(lengths(f) < 3)) die("a line in ", master, " has fewer than 3 tab-separated fields")
  d <- data.frame(slug = vapply(f, `[`, "", 1), doi = vapply(f, `[`, "", 2),
                  created_utc = vapply(f, `[`, "", 3),
                  title = vapply(f, function(z) if (length(z) >= 4) z[4] else "", ""),
                  stringsAsFactors = FALSE)
  if (any(!grepl(slug_ok, d$slug))) die("bad slug in ", master)
  if (any(!grepl(doi_any, d$doi)))  die("a DOI in ", master, " has an unexpected shape: ",
                                       paste(d$doi[!grepl(doi_any, d$doi)], collapse = ", "))
  if (anyDuplicated(d$slug)) die("duplicate slug in ", master, ": ", paste(d$slug[duplicated(d$slug)], collapse = ", "))
  if (anyDuplicated(d$doi))  die("duplicate DOI in ", master, ": ", paste(d$doi[duplicated(d$doi)], collapse = ", "))
  d
}
m <- read_master()

# ---- find the posts ---------------------------------------------------------
read_fm <- function(f) {
  x <- readLines(f, warn = FALSE, encoding = "UTF-8")
  if (!length(x) || trimws(x[1]) != "---") return(NULL)
  e <- which(trimws(x[-1]) %in% c("---", "..."))[1]
  if (is.na(e)) return(NULL)
  if (e < 2) return(list())
  tryCatch(yaml::yaml.load(paste(x[2:e], collapse = "\n")), error = function(err) NULL)
}
scalar <- function(v) if (is.null(v) || !length(v)) "" else paste(as.character(unlist(v)), collapse = " ")

idx   <- sort(Sys.glob("post/*/index.qmd"))
other <- setdiff(Sys.glob("post/*/*.qmd"), idx)
if (length(other)) warn("these .qmd files are not named index.qmd (they appear as separate listing items; ",
                        "this script does not handle them): ", paste(other, collapse = ", "))
cand <- data.frame(slug = basename(dirname(idx)), file = idx, title = "", status = "ok", stringsAsFactors = FALSE)
for (i in seq_len(nrow(cand))) {
  fm <- read_fm(cand$file[i])
  if (is.null(fm))                                  { cand$status[i] <- "skip: no readable front matter"; next }
  cand$title[i] <- gsub("[\t\r\n]+", " ", scalar(fm$title))
  if (tolower(scalar(fm$draft)) %in% c("true", "yes")) { cand$status[i] <- "skip: draft"; next }
  if (!grepl(slug_ok, cand$slug[i]))                { cand$status[i] <- "skip: slug has unusual characters"; next }
  if (!nzchar(scalar(fm$date)))                     { cand$status[i] <- "skip: no date (Rogue Scholar requires one)"; next }
  if (!nzchar(cand$title[i]))                       { cand$status[i] <- "skip: no title"; next }
}

todo    <- cand$slug[cand$status == "ok" & !(cand$slug %in% m$slug)]
skipped <- cand[cand$status != "ok", ]
orphans <- setdiff(m$slug, cand$slug)

info("Project: ", getwd())
info("Mode: ", if (write_mode) "WRITE (append to dois.tsv)" else "DRY RUN (nothing will be created)")
info("Posts found: ", nrow(cand), "   already have a DOI in ", master, ": ", sum(cand$slug %in% m$slug),
     "   need one: ", length(todo), "   skipped: ", nrow(skipped))
if (nrow(skipped)) for (i in seq_len(nrow(skipped))) info("  ", skipped$slug[i], "  [", skipped$status[i], "]")
if (length(orphans)) {
  warn(length(orphans), " slug(s) in ", master, " have no post folder now: ", paste(orphans, collapse = ", "))
  if (length(todo) && !allow_orphans)
    die("A slug with a DOI has disappeared while new slugs need DOIs.  That looks like a RENAMED post.\n",
        "Renaming a post folder changes its URL; do not do it (restore the old name).  If the old post was\n",
        "deliberately deleted, re-run with --allow-orphans.  Nothing was written.")
}
if (!length(todo)) { info("Nothing to do."); quit(save = "no", status = 0) }
if (!write_mode) {
  info("Would generate ", length(todo), " DOI(s) for:")
  for (s in todo) info("  ", s)
  info("Re-run with --write to do it.")
  quit(save = "no", status = 0)
}

# ---- generate and append, one line at a time --------------------------------
append_line <- function(line) {
  con <- file(master, open = "ab")
  on.exit(close(con))
  writeBin(charToRaw(enc2utf8(paste0(line, "\n"))), con)
  flush(con)
}
before <- if (file.exists(master)) readBin(master, "raw", n = file.size(master)) else raw(0)
if (!file.exists(master)) append_line(header)
used <- m$doi
n_ok <- 0L
for (s in sort(todo)) {
  d <- NULL
  for (try in 1:5) {
    out <- run_cm(c("encode", prefix)); st <- attr(out, "status")
    if (!is.null(st) && st != 0) die("commonmeta encode failed (status ", st, "); stopped after ", n_ok, " DOI(s)")
    cd <- norm_doi(tail(trimws(out[nzchar(trimws(out))]), 1))
    if (!grepl(doi_new, cd))  { warn("unexpected DOI shape from commonmeta: '", cd, "'"); next }
    if (!decode_ok(cd))       { warn("checksum validation failed for ", cd); next }
    if (cd %in% used)         { warn("duplicate of an existing DOI: ", cd); next }
    d <- cd; break
  }
  if (is.null(d)) die("could not get a valid, unique DOI for ", s, " after 5 tries; stopped after ", n_ok, " DOI(s)")
  tt <- cand$title[match(s, cand$slug)]
  append_line(paste(s, d, format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), tt, sep = "\t"))
  used <- c(used, d); n_ok <- n_ok + 1L
  info("  ", s, "  ", d)
}

# ---- self-check: the old content must be a byte-for-byte prefix of the new file
after <- readBin(master, "raw", n = file.size(master))
if (length(after) < length(before) || !identical(after[seq_along(before)], before))
  die("!!! the existing content of ", master, " changed.  Restore it from git NOW: git checkout -- ", master)
m2 <- read_master()
info("Appended ", n_ok, " DOI(s).  ", master, " now has ", nrow(m2), " entries.")
info("NEXT: git add ", master, " && git commit -m 'Add DOIs for ", n_ok, " post(s)'   (commit it immediately)")
