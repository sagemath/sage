# Cannot be found by pkg-config
symmetrica = cc.find_library('symmetrica')

py.install_sources('__init__.py', 'all.py', subdir: 'sage/libs/symmetrica')

extension_data = {'symmetrica' : files('symmetrica.pyx')}

foreach name, pyx : extension_data
  if is_windows
    # symmetrica/def.h(1769): error C4003: not enough arguments for function-like macro invocation 'max'
    # symmetrica/def.h(1769): error C2059: syntax error: '<parameter-list>'
    continue
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/symmetrica',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp, symmetrica],
  )
endforeach

