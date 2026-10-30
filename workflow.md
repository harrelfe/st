# fharrell.com (Statistical Thinking): DOI / Rogue Scholar workflow

Status as of 2026-10-29: 77 posts have DOIs (prefix 10.59350, registered via Rogue Scholar), active and resolving to https://www.fharrell.com/post/<slug>/. Deployed site verified with `check-live.sh` (ALL OK). Pending: ask Martin Fenner (Rogue Scholar) to refresh the 77 records after the feed changes below; then re-test the ecdf DOI in Zotero.

## Hard rules
- Never rename a slug. A DOI, once generated for a slug, never changes. `dois.tsv` (slug, doi, created_utc, title) is append-only; pre-commit hook `check-dois-append-only.sh` enforces it. Back up to ~/i/Backups/dois-YYYY-MM-DD.tsv.
- Canonical host is www.fharrell.com. Posts use `freeze: true` (do not re-run R code). Posts are CC BY 4.0.
- Every post must have its own complete `author:` list. The shared default was removed from `post/_metadata.yml` (Quarto concatenates author lists from directory metadata with the post's own, which wrongly added Frank Harrell to guest posts).

## Files (blog root, ~/i/R/blog/q)
- `make-dois.R` : dry run by default; `--write` appends new slugs to `dois.tsv` (uses commonmeta; self-tests checksum; needs stderr captured).
- `sync-doi-metadata.R` : `--write` adds `doi:` to `post/<slug>/_metadata.yml`.
- `cite-box.lua` (pre-quarto filter) : for posts with a DOI, puts a blue "How to cite this article" box at the top and adds head tags (citation_title, citation_author per author, citation_publication_date, citation_doi, DC.identifier, citation_abstract_html_url). Authors are read from the post's own .qmd front matter on disk (not the possibly-frozen copy); "A and B" in one name field is split; names may use given/family.
- `set-doi-guid.R` (post-render) : in `_site/index.xml` (a) drops non-/post/ items (talks), (b) removes the cite box from feed items, (c) inserts each post's front-matter `description:` as the first paragraph (`<p class="abstract-feed">`), (d) replaces `<guid>` with `https://doi.org/<doi>`. A full render stops (exit 1) if any post lacks a DOI. Idempotent. Summary line expected: removed 77 cite box(es); added 76 abstract(s); dropped 29 non-post item(s); replaced 77 post guid(s).
- `add-orcid.R`, `add-default-author.R` : one-time scripts, already used (ORCID added to the 72 posts Frank authored/co-authored).
- `check-live.sh` : run with `zsh check-live.sh` after deploy; checks live feed counts and, for every slug, page 200, canonical link, citation_doi, cite box, DOI guid in feed, and doi.org redirect target.
- `_quarto.yml` : `site-url`, `canonical-url: true` under `format: html:`, filter cite-box.lua, post-render `Rscript set-doi-guid.R` then syncdir.

## New-post workflow
1. Write the post with a full `author:` list (guests: use `name: {given: ..., family: ...}` for middle names or compound family names, else Quarto splits at the first space) and a `description:` (Rogue Scholar truncates abstracts near 500 chars, so put the key sentence first).
2. `Rscript make-dois.R` (check), then `--write`; `Rscript sync-doi-metadata.R --write`.
3. Full `quarto render` from the blog root (single-file renders re-run code and do not fully process the feed).
4. Check the summary line; commit; deploy; `zsh check-live.sh`.

## Facts learned
- Rogue Scholar abstracts come from the start of the feed item text; Quarto's full feed carries only the rendered body (not the front-matter description), hence the script's insertion step.
- Crossref registers posts as posted-content: blog title is in `group-title` ("Statistical Thinking"), publisher/institution is "Front Matter", container-title empty. Zotero's DOI import therefore shows Blog Title = "Front Matter"; edit by hand. Not fixable on our side.
- Corrections to registered records require emailing Martin (the feed has no updated date).
- Quarto `feed: false` on the talks listing is ignored; talks are removed by set-doi-guid.R.
- `google-scholar: false` in post/_metadata.yml (true crashed bibliography generation).
- macOS: BSD grep has no -P and caps `\{n,m\}` at 255; interactive zsh does not treat `#` as a comment; run scripts with `zsh`, not `sh`.
- Git/iCloud: repo lives in iCloud Drive (risk of duplicate "... 2" files; stale .git/index.lock fixed by moving it aside and `git reset`). `_site/post 2`, `talk 2`, `site_libs 2` duplicates should be cleaned.

## Open items
- Email Martin: refresh abstracts (cite box removed, descriptions added) and guest-post authors for the 77 records; ask what the earlier "minor issues with the content" were and how the DOI to archive redirect is triggered.
- Re-add the ecdf DOI in Zotero after the refresh.
- Fix the `ormll` post before editing it; consider moving the repo out of iCloud Drive.

## My Summary

cd ~/i/R/blog/q

_quarto.yml has  draft-mode: unlinked which means that draft articles are included in the generated html but are not findable by others who don't know the slug name of the article.

post/_metadata.yml has freeze:true so changes in non-R code will not run R

blog/q/index.qmd has feed: type: full  items: 1000; add feed: false for talk (also in talk/index.qmd)
blog/q/_quarto.yml has post-render: Rscript set-doi-guid.R and website: site-url: https://www.fharrell.com
Quarto would not respect feed: false or omitting feed, still producing feed for talks.  Modified set-doi-guid.R to remove talks from index.xml.

Have draft:true in article until it's really ready to publish/
To preview: open _site/index.html
Change URL to _site/post/slug/index.html
slug=directory name for article; article is in post/slug/index.*

When a new post is ready to be published, change draft to false and generate a permanent DOI for it:

1. wc -l dois.tsv
2. Rscript make-dois.R (dry run)
3. Rscript make-dois.R --write, which adds a line to dois.tsv.
4. wc -l dois.tsv
5. Rscript sync-doi-metadata.R (dry-run)
6. Rscript sync-doi-metadata.R --write (creates post/<slug>/_metadata.yml with doi: "...".)
7. cd slug and quarto render index.qmd or from blog root do quarto render
8. ./doi-backup
9. git add dois.tsv && git commit -m "add DOI" && gpush
10. _site_/index.xml to Martin of Rogue Scholar for initial checking
grep -c '<guid>https://doi.org/10.59350/' _site/index.xml
grep -c '<guid>' _site/index.xml
grep -c 'test-one' _site/index.xml
grep -c 'BY-SA' _site/index.xml

For each post listed in dois.tsv that has a post folder, sync-doi-metadata.R checks post/<slug>/_metadata.yml:

* File missing: it creates post/<slug>/_metadata.yml containing one line, doi: "10.59350/xxxxx-xxxxx".
* File exists with no doi: it appends that one line to the end, adding a newline first if the file didn’t end with one. Your other settings in the file are left alone.
* File already has the same DOI: it does nothing.
* File has a different DOI: it changes nothing and reports a conflict. dois.tsv is the master, so you fix the post file by hand.

The "how to cite" box at the top of each post is generated by 
1. filters: - at: pre-quarto path: cite-box.lua in _quarto.yml. This Lua filter runs on every page in the project
2. The metadata it reads. The filter builds the citation from each post’s front matter: title, author and date. It also uses the citation: block in post/_metadata.yml (type: webpage, container-title: "Statistical Thinking") and the post’s doi. The DOI comes from post/<slug>/_metadata.yml, which sync-doi-metadata.R --write fills in from dois.tsv.

The blog was added to the production Rogue Scholar site and the DOIs (which your blog already displays) are now registered: https://rogue-scholar.org/communities/fharrell

New blog posts are registered automatically. Send Martin Fenner martin@front-matter.de a short message when a post was updated and needs re-archiving.

Martin is making me owner of my blog community.  Posts are archived automatically, but you can change the community settings such as the logo.
