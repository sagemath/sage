doc_sources = [
  'advanced-2D-plotting.rst',
  'calculus.rst',
  'intro-tutorial.rst',
  'logging-on.rst',
  'programming.rst',
  'symbolics-and-basic-plotting.rst',
  'conf.py',
  'index.rst',
  'quickstart.rst',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('media')
subdir('quickstarts')
