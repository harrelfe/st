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


-- ---- authors -------------------------------------------------------------------------
local function trim(s) return (s:gsub("^%s+", ""):gsub("%s+$", "")) end

local function plain(x)
  if x == nil then return "" end
  if pandoc.utils.type(x) == "Inlines" then
    -- stop at the first raw HTML tag or line break (e.g. "Name<br><small>@handle</small>")
    local keep = pandoc.Inlines({})
    for _, el in ipairs(x) do
      if el.t == "RawInline" or el.t == "LineBreak" then break end
      keep:insert(el)
    end
    x = keep
  end
  return trim(pandoc.utils.stringify(x):gsub("<.*$", ""))   -- drop stray HTML such as a dangling "<br"
end

-- One author -> {family=, given=, full=}.  Handles "A B" strings, {name: "A B"},
-- {name: {given:, family:}} and {literal:}.
local parse_author
parse_author = function(a)
  local given, family
  local t = pandoc.utils.type(a)
  local nm = a
  if t == "table" or t == "Meta" or t == "Map" then
    nm = a.name or a.literal or a
  end
  if (pandoc.utils.type(nm) == "table" or pandoc.utils.type(nm) == "Map") and (nm.family or nm.given) then
    family, given = plain(nm.family), plain(nm.given)
  else
    local full = plain(nm)
    -- "A B and C D" in one name field means two people
    local pos = full:find(" and ", 1, true)
    if pos then
      local first  = parse_author(full:sub(1, pos - 1))
      local second = parse_author(full:sub(pos + 5))
      if first and second then
        for _, one in ipairs(second) do first[#first + 1] = one end
        return first
      end
    end
    local words = {}
    for w in full:gmatch("%S+") do words[#words + 1] = w end
    if #words == 0 then return nil end
    family = words[#words]
    given  = table.concat(words, " ", 1, #words - 1)
  end
  if family == "" then return nil end
  -- keep the middle initial the blog has always used for the owner
  if family == "Harrell" and (given == "Frank" or given == "Frank E." or given == "Frank E") then
    given = "Frank E."
  end
  return { { family = family, given = given } }
end

local function initials(given)
  local out = {}
  for w in given:gmatch("[^%s%.]+") do out[#out + 1] = pandoc.text.upper(pandoc.text.sub(w, 1, 1)) end
  return table.concat(out)
end

-- Authors come from the post's own .qmd front matter, read fresh from disk.  With `freeze: true`
-- the document pandoc sees is the frozen copy, whose YAML can be older than the .qmd.
local function fresh_author_meta()
  local infile = quarto and quarto.doc and quarto.doc.input_file or ""
  if infile == "" then return nil end
  local fh = io.open(infile, "r")
  if not fh then return nil end
  local lines, n, inyaml = {}, 0, false
  for line in fh:lines() do
    if line:match("^%-%-%-%s*$") then
      if not inyaml and n == 0 then inyaml = true
      else break end
    elseif inyaml then lines[#lines + 1] = line end
    n = n + 1
  end
  fh:close()
  if #lines == 0 then return nil end
  local ok, d = pcall(pandoc.read, "---\n" .. table.concat(lines, "\n") .. "\n---\n", "markdown")
  if ok and d and d.meta and d.meta.author then return d.meta end
  return nil
end

local function get_authors(m)
  local res = {}
  m = fresh_author_meta() or m
  local a = m.author
  if a == nil then return res end
  if pandoc.utils.type(a) ~= "List" then a = { a } end
  for _, x in ipairs(a) do
    local p = parse_author(x)
    if p then for _, one in ipairs(p) do res[#res + 1] = one end end
  end
  return res
end

function Pandoc(doc)
  local m = doc.meta
  if not m.doi then return nil end
  local doi   = pandoc.utils.stringify(m.doi)
  local title = pandoc.utils.stringify(m.title)
  local year, pubdate = parse_date(pandoc.utils.stringify(m.date or ""))
  local url   = "https://doi.org/" .. doi
  local authors = get_authors(m)
  if #authors == 0 then authors = { { family = "Harrell", given = "Frank E." } } end  -- safety net
  local cites = {}
  for _, au in ipairs(authors) do
    cites[#cites + 1] = au.family .. (au.given ~= "" and (" " .. initials(au.given)) or "")
  end
  local author_str = table.concat(cites, ", ")

  -- 1. visible box
  local heading = pandoc.Div(
    {pandoc.Plain({pandoc.Str("How to cite this article")})},
    pandoc.Attr("", {"cite-box-title"}))
  local para = pandoc.Para({
    pandoc.Str(author_str .. " (" .. year .. "). " .. title ..
               ". Statistical Thinking. "),
    pandoc.Link(url, url)
  })
  local box = pandoc.Div({heading, para}, pandoc.Attr("", {"cite-box"}))
  table.insert(doc.blocks, 1, box)

  -- 2. machine-readable tags in <head>
  local tags = {
    '<meta name="citation_title" content="' .. esc(title) .. '">',
    '<meta name="citation_publication_date" content="' .. esc(pubdate) .. '">',
    '<meta name="citation_doi" content="' .. esc(doi) .. '">',
    '<meta name="DC.identifier" content="' .. esc(url) .. '">',
  }
  for i, au in ipairs(authors) do
    table.insert(tags, 1 + i, '<meta name="citation_author" content="' ..
      esc(au.family .. (au.given ~= "" and (", " .. au.given) or "")) .. '">')
  end
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
