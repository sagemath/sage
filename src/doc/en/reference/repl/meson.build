doc_sources = [
  'conf.py',
  'environ.rst',
  'index.rst',
  'options.rst',
  'startup.rst',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

doc_src_repl = custom_target(
  'options',
  output: ['options.txt'],
  command: [py, src / 'sage' / 'cli', '--help'],
  capture: true,
)
