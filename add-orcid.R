#!/usr/bin/env Rscript
# add-orcid.R -- one-time fix.  Adds Frank Harrell's ORCID to every post whose front matter lists
# him as an author, because the shared `author:` block in post/_metadata.yml (the only place the
# ORCID currently lives) is being removed.
#
#   Rscript add-orcid.R           dry run
#   Rscript add-orcid.R --write   edit the files (use git to review)
#
# Matches an author line like   - name: Frank Harrell   (also "Frank E. Harrell", quoted or not)
# inside the front matter only, and skips authors that already have an orcid line.

write <- "--write" %in% commandArgs(TRUE)
orcid <- "0000-0002-8271-5493"
pat <- "^(\\s*)-\\s+name:\\s*[\"']?Frank (E\\.? )?Harrell[\"']?\\s*$"
n <- 0L
for (f in Sys.glob("post/*/index.qmd")) {
  l <- readLines(f, warn = FALSE, encoding = "UTF-8")
  e <- which(trimws(l) == "---")
  if (length(l) < 3 || e[1] != 1L || length(e) < 2) next
  i <- which(grepl(pat, l[1:e[2]]))
  if (!length(i)) next
  i <- i[1]
  if (any(grepl("^\\s*orcid:", l[(i + 1L):min(i + 6L, e[2])]))) next
  ind <- sub(pat, "\\1", l[i])
  n <- n + 1L
  cat(if (write) "ADD " else "would add ", f, "\n", sep = "")
  if (write) {
    l <- append(l, paste0(ind, "  orcid: ", orcid), after = i)
    con <- file(f, "wb"); writeLines(enc2utf8(l), con, useBytes = TRUE); close(con)
  }
}
cat(n, "post(s)", if (write) "changed" else "would change (dry run)", "\n")
