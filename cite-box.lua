function Pandoc(doc)
  local m = doc.meta
  if not m.doi then return nil end
  local doi   = pandoc.utils.stringify(m.doi)
  local title = pandoc.utils.stringify(m.title)
  local year  = pandoc.utils.stringify(m.date):match("%d%d%d%d") or ""
  local url   = "https://doi.org/" .. doi

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
  return doc
end
