-- Start these top-level sections on a fresh page (PDF output only).
local breaks = {
  Introduction = true, Methods = true, Results = true,
  Discussion = true, References = true,
}
function Header(h)
  if h.level == 1 and breaks[pandoc.utils.stringify(h)] then
    return { pandoc.RawBlock('latex', '\\clearpage'), h }
  end
end
