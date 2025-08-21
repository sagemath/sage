doc_src += custom_target(
  'doc-src',
  output: ['options.txt'],
  command: [py, src / 'sage' / 'cli', '--help'],
  capture: true,
)
