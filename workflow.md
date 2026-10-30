# fharrell.com (Statistical Thinking): DOI / Rogue Scholar workflow

**Status (2026-10-30).** 77 posts have permanent DOIs (prefix 10.59350). The blog is on the production Rogue Scholar site (<https://rogue-scholar.org/communities/fharrell>), the DOIs are registered and active, and each resolves to `https://www.fharrell.com/post/<slug>/`. The deployed site was verified with `check-live.sh` (ALL OK). New posts are registered automatically from the feed. Martin Fenner (martin@front-matter.de) is making Frank the owner of the Rogue Scholar community, which allows changing community settings such as the logo.

**Pending:** ask Martin to refresh the 77 existing records after the feed changes described below, then re-test the ecdf DOI in Zotero.

Project root: `~/i/R/blog/q` (all commands below run from there).

---

## 1. Publishing a new post: step by step

### 1.1 Draft and preview

- The article lives in `post/<slug>/index.qmd` (or another `index.*`). The **slug is the directory name**, and it never changes once published.
- Keep `draft: true` in the article until it is really ready to publish. `_quarto.yml` has `draft-mode: unlinked`, so drafts are rendered into the site but are not linked anywhere; only someone who knows the slug can find them.
- Give the post a **full `author:` list** (there is no shared default) and a **`description:`** (see section 4).
- To preview: open `_site/index.html`, or change the URL to `_site/post/<slug>/index.html`.
- `post/_metadata.yml` has `freeze: true`, so changes to non-R text do not re-run R code when the whole project is rendered.

### 1.2 Publish: generate the permanent DOI

When the post is ready, change `draft: true` to `draft: false`, then:

1. `wc -l dois.tsv`
2. `Rscript make-dois.R` (dry run)
3. `Rscript make-dois.R --write` (adds one line to `dois.tsv`)
4. `wc -l dois.tsv` (should be one more than in step 1)
5. `Rscript sync-doi-metadata.R` (dry run)
6. `Rscript sync-doi-metadata.R --write` (creates or updates `post/<slug>/_metadata.yml` with `doi: "..."`)
7. **Render the whole project from the blog root:** `quarto render`. Check the summary line printed by `set-doi-guid.R`, for example `removed 77 cite box(es); added 76 abstract(s); dropped 29 non-post item(s); replaced 77 post guid(s)` (the numbers grow as posts are added). A full render stops with an error if any post in the feed has no DOI.
   - Rendering only the post (`cd post/<slug>; quarto render index.qmd`) is fine for previewing, but it re-runs the R code and does not fully process the feed, so it is not enough before publishing.
8. `./doi-backup` (copies `dois.tsv` to `~/i/Backups/dois-YYYY-MM-DD.tsv`)
9. `git add dois.tsv && git commit -m "add DOI" && gpush`. Also commit the new post folder, including `post/<slug>/_metadata.yml`, either in the same commit or in another.
10. Deploy as usual.
11. `zsh check-live.sh` (run with `zsh`, not `sh`). It checks the live feed and, for every slug in `dois.tsv`, that the page loads, has the canonical link, `citation_doi` tag and cite box, that the feed carries the DOI as the item's guid, and that `doi.org/<doi>` redirects to the post. Expect `ALL OK`; an occasional HTTP `000` is a dropped connection, so rerun.

### 1.3 Feed sanity checks (optional, after step 7)

```
grep -c '<guid>https://doi.org/10.59350/' _site/index.xml
grep -c '<guid>' _site/index.xml
grep -c 'test-one' _site/index.xml
grep -c 'BY-SA' _site/index.xml
```

The first two counts should be equal (the number of posts). The last two should be `0` (no test post, and the license must be CC BY, not BY-SA).

### 1.4 After the post is published

- New posts are registered by Rogue Scholar automatically from the feed.
- If a published post is **updated** and needs re-archiving, send Martin a short message. The feed carries no updated date, so he cannot detect the change himself. The same applies to corrections of registered metadata.

---

## 2. How the pieces fit

### 2.1 Feed and site configuration

- `_quarto.yml`: `website: site-url: https://www.fharrell.com`; `canonical-url: true` (under `format: html:`, not `website:`); the filter `- at: pre-quarto, path: cite-box.lua`; and `post-render: Rscript set-doi-guid.R`, followed by `syncdir _site/ /Users/harrelfe/web/blog`.
- `index.qmd` (blog root): the post listing has `feed: type: full, items: 1000`. `feed: false` was also added for the talk listing (and in `talk/index.qmd`), but Quarto ignores it and still puts talks in the feed. `set-doi-guid.R` therefore removes talks from `index.xml`.
- `post/_metadata.yml`: `freeze: true`, `google-scholar: false` (`true` crashed bibliography generation), license CC BY, and the `citation:` block (`type: webpage`, `container-title: "Statistical Thinking"`). There is deliberately **no** `author:` here.

### 2.2 Scripts (blog root)

| File | Purpose |
|---|---|
| `make-dois.R` | Dry run by default; `--write` appends new slugs to `dois.tsv`. Uses commonmeta, self-tests the DOI checksum. |
| `sync-doi-metadata.R` | Copies each DOI from `dois.tsv` into `post/<slug>/_metadata.yml` (rules in 2.3). |
| `check-dois-append-only.sh` | Git pre-commit hook: `dois.tsv` may only grow. |
| `cite-box.lua` | Pre-quarto filter: the "How to cite this article" box and the citation head tags (2.4). |
| `set-doi-guid.R` | Post-render feed fixer (2.5). |
| `check-live.sh` | Post-deploy verification (step 11). |
| `doi-backup` | Backs up `dois.tsv` to `~/i/Backups`. |
| `add-orcid.R`, `add-default-author.R` | One-time scripts, already used. ORCID was added to the 72 posts Frank wrote or co-wrote. |

### 2.3 What `sync-doi-metadata.R` does

For each post listed in `dois.tsv` that has a post folder, it checks `post/<slug>/_metadata.yml`:

- **File missing:** creates it containing one line, `doi: "10.59350/xxxxx-xxxxx"`.
- **File exists with no `doi`:** appends that one line (adding a newline first if the file didn't end with one). Other settings in the file are left alone.
- **File already has the same DOI:** does nothing.
- **File has a different DOI:** changes nothing and reports a conflict. `dois.tsv` is the master, so fix the post file by hand.

### 2.4 The "How to cite" box and citation tags (`cite-box.lua`)

The filter (`filters: - at: pre-quarto, path: cite-box.lua` in `_quarto.yml`) runs on every page in the project but only acts on pages that have a `doi`. The DOI comes from `post/<slug>/_metadata.yml`, which `sync-doi-metadata.R --write` fills in from `dois.tsv`.

- It builds the citation line from the post's title, date and authors, plus the DOI link, in the form "Levy DG, Steyerberg E (2018). Title. Statistical Thinking. https://doi.org/...". The blog name "Statistical Thinking" is fixed in the filter. (The `citation:` block in `post/_metadata.yml` is used by Quarto itself for the BibTeX and attribution at the bottom of the page, not by this filter.)
- It puts the blue box (CSS in `h.css`) at the top of the page and adds head tags: `citation_title`, one `citation_author` per author, `citation_publication_date`, `citation_doi`, `DC.identifier` and `citation_abstract_html_url`.
- **Authors are read from the post's own `.qmd` file on disk**, not from the pandoc copy, which may be an older frozen version. A name field containing "A and B" is split into two authors.

### 2.5 What `set-doi-guid.R` does to `_site/index.xml`

1. Drops every item whose link is not under `/post/` (talks).
2. Removes the cite box from each item. It is useful on the web page, but Rogue Scholar builds each abstract from the start of the item text, so the box would become the start of every abstract.
3. Inserts each post's front-matter `description:` as the first paragraph of its item (`<p class="abstract-feed">`). Quarto's full-text feed carries only the rendered body, not the description.
4. Replaces each `<guid>` (the post URL) with `https://doi.org/<doi>`. This DOI URL is what Rogue Scholar uses as the post's permanent identifier.

A full render stops with exit 1 if any post has no DOI. The script is idempotent: running it twice changes nothing.

---

## 3. Hard rules

- **Never rename a slug.** A DOI, once generated for a slug, never changes.
- **`dois.tsv` is append-only** (columns: slug, doi, created_utc, title). The pre-commit hook enforces it. Back it up after every addition.
- **The canonical host is `www.fharrell.com`.**
- **Never re-run R code in posts** (`freeze: true`).
- Posts are licensed **CC BY 4.0** (Rogue Scholar requires CC BY).
- **Every post has its own complete `author:` list.** Quarto concatenates the `author:` list from `post/_metadata.yml` with the post's own list, which had added Frank Harrell to guest posts, so the shared default was removed.
  - Guest authors with a middle name or compound family name must be written `name: {given: Drew Griffin, family: Levy}`. Otherwise Quarto splits at the first space and puts "Griffin Levy" in the family name.
  - Write each author as a separate entry; do not use "A and B" in one `name:`.

---

## 4. What Rogue Scholar and Crossref do with the posts

- Rogue Scholar reads the RSS feed (`_site/index.xml`), archives the posts, and registers the DOIs with Crossref. Corrections to registered records go through Martin.
- **Abstracts** come from the start of the feed item text, truncated near 500 characters. Put the key sentence first in each `description:`. Updating the abstract of an existing record requires Martin to re-import.
- Crossref registers the posts as "posted content": the blog title ("Statistical Thinking") is in `group-title`, the publisher and institution are "Front Matter", and `container-title` is empty. **Zotero's DOI import therefore shows Blog Title = "Front Matter".** Edit it by hand. It cannot be fixed on our side.
- A post that is updated after publication keeps its DOI; tell Martin so it can be re-archived (1.4).

---

## 5. Gotchas

- **macOS:** BSD `grep` has no `-P` and limits `\{n,m\}` repetition to 255. Interactive zsh does not treat `#` as a comment on the command line. Run scripts with `zsh`, not `sh` (`sh` lacks `< <(...)`).
- **Quarto:** `canonical-url` belongs under `format: html:` and needs `site-url`. Only a full render from the blog root (not a single-file or subdirectory render) fully processes the feed.
- **Git and iCloud:** the repository lives in iCloud Drive, which risks duplicate files named "... 2" (for example `_site/post 2`, `talk 2`, `site_libs 2`; these should be cleaned up). A stale `.git/index.lock` makes git show every file as a staged deletion plus untracked; fix it by moving the lock aside and running `git reset`.
- `_book` in another repo (the R Workflow book): tracked files ignore `.gitignore`; use `git rm -r --cached` and check with `git check-ignore --no-index -v`.

---

## 6. Open items

- Email Martin: refresh the 77 records' abstracts (cite box removed, descriptions added) and the author lists of the guest posts. Ask what the earlier "minor issues with the content" were, and how the DOI-to-archive redirect is triggered.
- Re-add the ecdf DOI in Zotero after the refresh and check the abstract and authors.
- Fix the `ormll` post before editing it.
- Consider moving the repository out of iCloud Drive; clean up the "... 2" duplicates.
