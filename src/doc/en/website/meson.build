doc_sources = ['conf.py', 'index.rst', 'root_index.html', 'versions.txt']

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('templates')
