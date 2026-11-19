local function has_class(classes, name)
  for _, class in ipairs(classes) do
    if class == name then
      return true
    end
  end
  return false
end

local function render_figure(src, caption, fullpage)
  local height = fullpage and '0.88\\textheight' or '0.52\\textheight'
  local blocks = {}

  if fullpage then
    blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\clearpage')
    blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\newgeometry{margin=0.5in}')
    src = src:gsub('%.png$', '.pdf')
  end

  blocks[#blocks + 1] = pandoc.RawBlock('latex',
    '\\begin{center}\n' ..
    '\\includegraphics[' ..
    'width=\\linewidth,' ..
    'height=' .. height .. ',' ..
    'keepaspectratio' ..
    ']{\\detokenize{' .. src .. '}}\n' ..
    '\\end{center}')

  if #caption > 0 then
    blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\begingroup\\small')
    for _, block in ipairs(caption) do
      blocks[#blocks + 1] = block
    end
    blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\endgroup')
  end

  if fullpage then
    blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\clearpage')
    blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\restoregeometry')
  end

  return blocks
end

function Meta(meta)
  if FORMAT:match('latex') then
    local includes = {}
    local existing = meta['header-includes']
    if existing then
      if existing.t == 'MetaBlocks' then
        includes[#includes + 1] = existing
      else
        for _, include in ipairs(existing) do
          includes[#includes + 1] = include
        end
      end
    end
    includes[#includes + 1] = pandoc.MetaBlocks({
      pandoc.RawBlock('latex', '\\usepackage{graphicx}\\usepackage{geometry}')
    })
    meta['header-includes'] = pandoc.MetaList(includes)
  end
  return meta
end

function Figure(el)
  if not FORMAT:match('latex') then
    return nil
  end

  for _, block in ipairs(el.content) do
    if block.content then
      for _, inline in ipairs(block.content) do
        if inline.t == 'Image' then
          local fullpage = has_class(el.attr.classes, 'fullpage') or
            has_class(inline.classes, 'fullpage')
          return render_figure(inline.src, el.caption.long or {}, fullpage)
        end
      end
    end
  end

  return nil
end

function Para(el)
  if FORMAT:match('latex') then
    local caption = pandoc.utils.stringify(el)
    if caption:match('^Table%s+%d+%.') then
      table.insert(el.content, 1, pandoc.RawInline('latex', '\\begingroup\\small '))
      table.insert(el.content, pandoc.RawInline('latex', '\\endgroup'))
      return el
    end
  end
end
