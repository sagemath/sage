bliss = cc.find_library('bliss', required: get_option('extension_bliss_standalone'), disabler: true)

py.extension_module(
  '__init__',
  sources: files('__init__.pyx'),
  subdir: 'sage/graphs/bliss',
  install: true,
  override_options: ['cython_language=cpp'],
  include_directories: [],
  dependencies: [cysignals, bliss],
)
