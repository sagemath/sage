symmetrica = dependency(
  'symmetrica',
  required: get_option('symmetrica'),
  disabler: true,
)
conf_data.set10('SYMMETRICA_ENABLED', symmetrica.found())

py.install_sources(
  '__init__.py',
  'all.py',
  'symmetrica.pyx',
  subdir: 'sage/libs/symmetrica',
)

extension_data = {'symmetrica': files('symmetrica.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/symmetrica',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp, symmetrica],
  )
endforeach

