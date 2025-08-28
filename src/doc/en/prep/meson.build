doc_sources = [
  'Advanced-2DPlotting.rst',
  'Calculus.rst',
  'Intro-Tutorial.rst',
  'Logging-On.rst',
  'Programming.rst',
  'Symbolics-and-Basic-Plotting.rst',
  'conf.py',
  'index.rst',
  'quickstart.rst',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('media')
subdir('Quickstarts')
