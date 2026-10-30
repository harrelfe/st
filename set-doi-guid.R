#!/usr/bin/env Rscript
# set-doi-guid.R -- Quarto post-render step.  In the generated RSS feed, replace each
# post's <guid> (Quarto uses the post URL) with the post's DOI URL, taken from the
# master file dois.tsv.  Rogue Scholar uses the guid as the post's permanent identifier.
#
# Quarto supplies:  QUARTO_PROJECT_OUTPUT_DIR (default _site)
#                   QUARTO_PROJECT_RENDER_ALL (set on a full-project render)
# Optional:         DOI_DRYRUN=1  leave the feed alone; write index-doi-preview.xml in the project root
#
# Behavior
#   * Only guids of the form  <guid>https://HOST/post/SLUG/</guid>  are touched.  Talks, and
#     guids already converted to a DOI, are left alone, so running it twice changes nothing.
#   * Full render: if any post in the feed has no DOI in dois.tsv, the script STOPS (exit 1)
#     before writing anything, so a post can never be published to Rogue Scholar with a random DOI.
#     Single-file render: only warns.
#   * Items whose <link> is not under /post/ (talks, etc.) are removed from the feed.
#   * Every other byte of the feed outside the replaced guid text is preserved.

die  <- function(...) { cat("ERROR:", ..., "\n", file = stderr()); quit(save = "no", status = 1) }
info <- function(...) cat(..., "\n", sep = "")

outdir <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR", "_site")
feed   <- file.path(outdir, "index.xml")
full   <- nzchar(Sys.getenv("QUARTO_PROJECT_RENDER_ALL"))
dry    <- nzchar(Sys.getenv("DOI_DRYRUN"))
master <- "dois.tsv"

if (!file.exists(feed)) quit(save = "no", status = 0)        # e.g. a single-post render: no feed to fix
if (!file.exists(master)) die(master, " not found in ", getwd())

ln <- readLines(master, warn = FALSE, encoding = "UTF-8")
if (!length(ln) || ln[1] != "slug\tdoi\tcreated_utc\ttitle") die("unexpected header in ", master)
f <- strsplit(ln[-1], "\t", fixed = TRUE)
if (any(lengths(f) < 3)) die("malformed line in ", master)
slug <- vapply(f, `[`, "", 1)
doi  <- tolower(sub("^https?://(dx\\.)?doi\\.org/", "", vapply(f, `[`, "", 2)))
if (anyDuplicated(slug) || anyDuplicated(doi)) die("duplicate slug or DOI in ", master)
if (any(!grepl("^10\\.[0-9]{4,9}/[0-9a-z]{5}-[0-9a-z]{3}[0-9]{2}$", doi))) die("malformed DOI in ", master)
map <- setNames(doi, slug)

raw <- readBin(feed, "raw", n = file.size(feed))
x   <- rawToChar(raw); Encoding(x) <- "UTF-8"

# Drop every feed item whose <link> is not under /post/ (e.g. talks): Rogue Scholar registers
# everything in the feed permanently, and only blog posts belong there.  Items without a
# <link> are kept.  Bytes outside the dropped items are untouched.
st <- gregexpr("<item>",  x, fixed = TRUE, useBytes = TRUE)[[1]]
en <- gregexpr("</item>", x, fixed = TRUE, useBytes = TRUE)[[1]]
n_dropped <- 0L
if (st[1] > 0 && length(st) == length(en) && all(en > st) && all(st[-1] > en[-length(en)])) {
  keep <- rep(TRUE, length(raw))
  for (i in seq_along(st)) {
    it <- rawToChar(raw[st[i]:(en[i] + 6L)]); Encoding(it) <- "UTF-8"
    lk <- regmatches(it, regexpr("<link>[^<]*</link>", it, useBytes = TRUE))
    if (length(lk) && !grepl("/post/", lk, fixed = TRUE)) {
      a <- st[i]
      while (a > 1 && raw[a - 1] %in% as.raw(c(0x20, 0x09))) a <- a - 1L      # indentation before <item>
      b <- en[i] + 6L
      if (b < length(raw) && raw[b + 1] == as.raw(0x0a)) b <- b + 1L           # newline after </item>
      keep[a:b] <- FALSE; n_dropped <- n_dropped + 1L
    }
  }
  if (n_dropped) {
    raw <- raw[keep]
    x   <- rawToChar(raw); Encoding(x) <- "UTF-8"
  }
} else if (st[1] > 0) warning("feed <item> tags are not well paired; non-post items were not filtered", call. = FALSE)
# Remove the "How to cite this article" box (added by cite-box.lua) from the feed items.  It is
# useful on the web page, but Rogue Scholar builds each record's abstract from the start of the
# item text, so leaving it in puts the citation into every abstract.  The pages are not touched.
pat_cb <- "(?s)<div class=\"cite-box\">\\s*<div class=\"cite-box-title\">.*?</div>\\s*<p>.*?</p>\\s*</div>\\s*"
n_cite <- length(gregexpr(pat_cb, x, perl = TRUE, useBytes = TRUE)[[1]])
if (gregexpr(pat_cb, x, perl = TRUE, useBytes = TRUE)[[1]][1] < 0) n_cite <- 0L
x <- gsub(pat_cb, "", x, perl = TRUE, useBytes = TRUE)
# Put each post's front-matter `description:` at the start of its feed item.  Quarto's full-text
# feed carries only the rendered body, so without this Rogue Scholar builds the abstract from
# whatever the post happens to start with.  Posts with no description are left alone, and items
# that already have the abstract paragraph are skipped, so a second run changes nothing.
Encoding(x) <- "UTF-8"
get_desc <- function(slug) {
  q <- file.path("post", slug, "index.qmd")
  if (!file.exists(q)) return("")
  l <- readLines(q, warn = FALSE, encoding = "UTF-8")
  e <- which(trimws(l) == "---")
  if (length(e) < 2 || e[1] != 1L) return("")
  d <- tryCatch(yaml::yaml.load(paste(l[(e[1] + 1L):(e[2] - 1L)], collapse = "\n"))$description,
                error = function(err) NULL)
  if (!is.character(d) || length(d) != 1L) return("")
  d <- gsub("---", "\u2014", d, fixed = TRUE); d <- gsub("--", "\u2013", d, fixed = TRUE)
  d <- gsub("&", "&amp;", d, fixed = TRUE); d <- gsub("<", "&lt;", d, fixed = TRUE); d <- gsub(">", "&gt;", d, fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", d))
}
n_abs <- 0L
its <- regmatches(x, gregexpr("(?s)<item>.*?</item>", x, perl = TRUE))[[1]]
if (length(its)) {
  new <- vapply(its, function(it) {
    sl <- sub("(?s).*<link>[^<]*/post/([^/<]+)/?</link>.*", "\\1", it, perl = TRUE)
    if (identical(sl, it) || grepl('class="abstract-feed"', it, fixed = TRUE)) return(it)
    d <- get_desc(sl)
    if (!nzchar(d)) return(it)
    n_abs <<- n_abs + 1L
    sub("(<description><!\\[CDATA\\[\\s*)", paste0("\\1<p class=\"abstract-feed\">", gsub("\\", "\\\\", d, fixed = TRUE), "</p>\n"), it, perl = TRUE)
  }, "", USE.NAMES = FALSE)
  regmatches(x, gregexpr("(?s)<item>.*?</item>", x, perl = TRUE)) <- list(new)
}
n_doi_before <- length(regmatches(x, gregexpr("<guid>https://doi.org/", x, fixed = TRUE, useBytes = TRUE))[[1]])

pat <-"<guid>[^<]*/post/[^/<]+/?</guid>"
g   <- unique(regmatches(x, gregexpr(pat, x, useBytes = TRUE))[[1]])
s   <- sub("^<guid>.*/post/([^/<]+)/?</guid>$", "\\1", g)
have <- s %in% names(map)

if (any(!have)) {
  msg <- paste0(sum(!have), " post(s) in the feed have no DOI in ", master, ": ", paste(s[!have], collapse = ", "))
  if (full) die(msg, "\nNothing was changed.  Run: Rscript make-dois.R --write")
  warning(msg, call. = FALSE)
}
for (i in which(have)) {
  x <- sub(g[i], sprintf("<guid>https://doi.org/%s</guid>", map[[s[i]]]), x, fixed = TRUE)
}

out <- if (dry) "index-doi-preview.xml" else feed
writeBin(charToRaw(x), out)
info("guid -> DOI: removed ", n_cite, " cite box(es); added ", n_abs, " abstract(s); dropped ", n_dropped, " non-post item(s); replaced ", sum(have), " post guid(s); ", n_doi_before, " were already DOIs; ",
     length(regmatches(x, gregexpr("<guid>", x, fixed = TRUE, useBytes = TRUE))[[1]]), " guid(s) total in the feed",
     if (dry) paste0("  [DRY RUN: wrote ", out, "]") else "")
