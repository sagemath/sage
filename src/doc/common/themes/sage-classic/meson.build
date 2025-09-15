doc_sources = [
  'genindex-single.html',
  'genindex-split.html',
  'genindex.html',
  'layout.html',
  'search.html',
  'theme.conf',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('static')
