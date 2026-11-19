local function render_figure(src, caption)
  local blocks = {}

  blocks[#blocks + 1] = pandoc.RawBlock('latex',
    '\\begin{center}\n' ..
    '\\includegraphics[' ..
    'width=\\linewidth,' ..
    'height=0.62\\textheight,' ..
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

  return blocks
end

local function table_widths(table)
  local column_count = #table.colspecs
  local presets = {
    [3] = {0.40, 0.30, 0.30},
    [5] = {0.34, 0.16, 0.16, 0.16, 0.18},
    [6] = {0.13, 0.28, 0.1475, 0.1475, 0.1475, 0.1475},
    [7] = {0.28, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12},
    [8] = {0.13, 0.17, 0.14, 0.10, 0.14, 0.14, 0.08, 0.10},
  }
  if column_count == 7 and tostring(table.colspecs[1][1]) == 'AlignRight' then
    return {0.05, 0.10, 0.16, 0.43, 0.08, 0.08, 0.10}
  end
  return presets[column_count]
end

local function size_table(table)
  local widths = table_widths(table)
  if widths then
    for index, colspec in ipairs(table.colspecs) do
      table.colspecs[index] = {colspec[1], widths[index]}
    end
  end
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
      pandoc.RawBlock('latex',
        '\\usepackage{graphicx}\\usepackage{geometry}' ..
        '\\usepackage{pdflscape}\\usepackage{needspace}')
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
          return render_figure(inline.src, el.caption.long or {})
        end
      end
    end
  end

  return nil
end

function Pandoc(doc)
  if not FORMAT:match('latex') then
    return doc
  end

  local blocks = {}
  local index = 1
  while index <= #doc.blocks do
    local block = doc.blocks[index]
    local next_block = doc.blocks[index + 1]
    local is_table_caption = block.t == 'Para' and
      pandoc.utils.stringify(block):match('^Table%s+%d+%.')

    if is_table_caption and next_block and next_block.t == 'Table' then
      local columns = #next_block.colspecs
      size_table(next_block)
      local needed_space = columns == 6 and 20 or 4
      blocks[#blocks + 1] = pandoc.RawBlock('latex',
        '\\needspace{' .. needed_space .. '\\baselineskip}')

      if columns >= 7 then
        blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\begin{landscape}')
      end

      blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\begingroup\\small')
      blocks[#blocks + 1] = block
      blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\endgroup')
      blocks[#blocks + 1] = next_block

      if columns >= 7 then
        blocks[#blocks + 1] = pandoc.RawBlock('latex', '\\end{landscape}')
      end
      index = index + 2
    else
      blocks[#blocks + 1] = block
      index = index + 1
    end
  end

  doc.blocks = blocks
  return doc
end
