doc_sources = ['conf.py', 'index.rst', 'threejs.rst']

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('threejs_examples')
