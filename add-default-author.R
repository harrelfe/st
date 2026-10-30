#!/usr/bin/env Rscript
# Note: this is not needed since all blog posts already had author:
# add-default-author.R -- one-time fix.  Quarto concatenates the `author:` list in
# post/_metadata.yml with each post's own `author:` list, so posts by other people also list
# Frank Harrell.  After you delete the `author:` block from post/_metadata.yml, every post that
# has no author of its own would have none at all; this adds Frank Harrell to exactly those.
#
#   Rscript add-default-author.R           dry run: list the files that would change
#   Rscript add-default-author.R --write   edit them (a .bak copy is not made: use git)
#
# A post is left alone if its front matter already has `author:` or `authors:`.

write <- "--write" %in% commandArgs(TRUE)
block <- c("author:",
           "  - name: Frank Harrell",
           "    orcid: 0000-0002-8271-5493",
           "    affiliation: Vanderbilt University School of Medicine")
files <- Sys.glob("post/*/index.qmd")
n_chg <- 0L
for (f in files) {
  l <- readLines(f, warn = FALSE, encoding = "UTF-8")
  e <- which(trimws(l) == "---")
  if (length(l) < 3 || e[1] != 1L || length(e) < 2) { cat("SKIP (no front matter):", f, "\n"); next }
  fm <- l[2:(e[2] - 1L)]
  if (any(grepl("^(author|authors):", fm))) next
  n_chg <- n_chg + 1L
  cat(if (write) "ADD " else "would add ", f, "\n", sep = "")
  if (write) {
    l <- c(l[1], block, l[-1])
    con <- file(f, "wb"); writeLines(enc2utf8(l), con, useBytes = TRUE); close(con)
  }
}
cat(n_chg, "post(s)", if (write) "changed" else "would change (dry run)", "out of", length(files), "\n")
