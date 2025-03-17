# Cannot be found by pkg-config
symmetrica = cc.find_library('symmetrica')

py.install_sources('__init__.py', 'all.py', subdir: 'sage/libs/symmetrica')

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

