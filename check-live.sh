#!/bin/zsh
# check-live.sh -- verify the deployed blog against dois.tsv.  Run from the blog root.
# For every slug: the page loads, has the canonical link, citation_doi and cite box; the DOI
# redirects to the page; the live feed carries the DOI as the item's guid.
base=https://www.fharrell.com
tmp=$(mktemp -d)
fail=0
c=(curl --retry 4 --retry-delay 2 --retry-all-errors --connect-timeout 15 --max-time 60)

$c -sL -o $tmp/feed.xml $base/index.xml || { echo "cannot fetch feed"; exit 1; }
echo "== live feed =="
echo "cite-box occurrences (want 0):       $(grep -c 'cite-box' $tmp/feed.xml)"
echo "abstract-feed paragraphs:            $(grep -c 'abstract-feed' $tmp/feed.xml)"
echo "DOI guids (want number of posts):    $(grep -c '<guid>https://doi.org/' $tmp/feed.xml)"
echo "post-URL guids (want 0):             $(grep -c '<guid>https://www.fharrell.com/post/' $tmp/feed.xml)"
echo "talk items (want 0):                 $(grep -c '<link>[^<]*/talk/' $tmp/feed.xml)"
if command -v xmllint >/dev/null; then
  xmllint --noout $tmp/feed.xml && echo "feed is well-formed XML" || { echo "FEED NOT WELL-FORMED"; fail=1; }
fi

echo "\n== per post (only problems are printed) =="
n=0
while IFS=$'\t' read -r slug doi rest; do
  [[ -z "$slug" ]] && continue
  n=$((n+1))
  doi=${doi#https://doi.org/}
  url=$base/post/$slug/
  code=$($c -sL -o $tmp/p.html -w '%{http_code}' $url)
  [[ $code == 200 ]] || { echo "$slug: page returned HTTP $code (a code of 000 is a connection failure, not a site problem; rerun)"; fail=1; continue; }
  grep 'rel="canonical"' $tmp/p.html | grep -q "$url"        || { echo "$slug: canonical link missing or wrong"; fail=1; }
  grep -q "name=\"citation_doi\" content=\"$doi\"" $tmp/p.html || { echo "$slug: citation_doi missing or wrong"; fail=1; }
  grep -q 'class="cite-box"' $tmp/p.html                       || { echo "$slug: cite box missing"; fail=1; }
  grep -q "<guid>https://doi.org/$doi</guid>" $tmp/feed.xml    || { echo "$slug: DOI guid not in feed"; fail=1; }
  red=$($c -s -o /dev/null -w '%{redirect_url}' https://doi.org/$doi)
  [[ $red == $url ]] || { echo "$slug: DOI redirects to '$red' (want $url)"; fail=1; }
done < <(tail -n +2 dois.tsv)

echo "\nchecked $n posts; " $([[ $fail == 0 ]] && echo "ALL OK" || echo "PROBLEMS FOUND")
rm -rf $tmp
