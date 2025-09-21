symmetrica = dependency('symmetrica', required: false, disabler: true)
if not symmetrica.found()
  # fallback for older versions without pkg-config
  symmetrica = cc.find_library(
    'symmetrica',
    required: not is_windows,
    disabler: true,
  )
endif

py.install_sources(
  '__init__.py',
  'all.py',
  'symmetrica.pyx',
  subdir: 'sage/libs/symmetrica',
)

extension_data = {'symmetrica' : files('symmetrica.pyx')}

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

