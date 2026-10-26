-- cite-box.lua: for posts that have a `doi`, (1) put a "How to cite this article" box at the top
-- and (2) add citation / Dublin Core meta tags to the page <head>.

local months = { january=1, february=2, march=3, april=4, may=5, june=6, july=7,
                 august=8, september=9, october=10, november=11, december=12 }

local function esc(s)
  return (s:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;"):gsub('"', "&quot;"))
end

-- Returns year, "YYYY/MM/DD" (or just year if the day/month can't be read)
local function parse_date(d)
  local y, mo, dd = d:match("(%d%d%d%d)[-/](%d%d?)[-/](%d%d?)")
  if y then return y, string.format("%s/%02d/%02d", y, tonumber(mo), tonumber(dd)) end
  local mname, day, yr = d:match("(%a+)%s+(%d%d?),?%s+(%d%d%d%d)")
  if mname and months[mname:lower()] then
    return yr, string.format("%s/%02d/%02d", yr, months[mname:lower()], tonumber(day))
  end
  local year = d:match("%d%d%d%d") or ""
  return year, year
end

function Pandoc(doc)
  local m = doc.meta
  if not m.doi then return nil end
  local doi   = pandoc.utils.stringify(m.doi)
  local title = (pandoc.utils.stringify(m.title):gsub("\u{00A0}", " "))
  local year, pubdate = parse_date(pandoc.utils.stringify(m.date or ""))
  local url   = "https://doi.org/" .. doi

  -- 1. visible box
  local heading = pandoc.Div(
    {pandoc.Plain({pandoc.Str("How to cite this article")})},
    pandoc.Attr("", {"cite-box-title"}))
  local para = pandoc.Para({
    pandoc.Str("Harrell FE (" .. year .. "). " .. title ..
               ". Statistical Thinking. "),
    pandoc.Link(url, url)
  })
  local box = pandoc.Div({heading, para}, pandoc.Attr("", {"cite-box"}))
  table.insert(doc.blocks, 1, box)

  -- 2. machine-readable tags in <head>
  local tags = {
    '<meta name="citation_title" content="' .. esc(title) .. '">',
    '<meta name="citation_author" content="Harrell, Frank E.">',
    '<meta name="citation_publication_date" content="' .. esc(pubdate) .. '">',
    '<meta name="citation_doi" content="' .. esc(doi) .. '">',
    '<meta name="DC.identifier" content="' .. esc(url) .. '">',
  }
  local infile = quarto and quarto.doc and quarto.doc.input_file or ""
  local slug = infile:match("post/([^/]+)/index%.qmd$")
  if slug then
    tags[#tags + 1] = '<meta name="citation_abstract_html_url" content="https://www.fharrell.com/post/'
                      .. esc(slug) .. '/">'
  end
  if quarto and quarto.doc and quarto.doc.include_text then
    quarto.doc.include_text("in-header", table.concat(tags, "\n"))
  end

  return doc
end
