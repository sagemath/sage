# Cannot be found by pkg-config
symmetrica = cc.find_library('symmetrica')

py.install_sources('__init__.py', 'all.py', subdir: 'sage/libs/symmetrica')

extension_data = {'symmetrica' : files('symmetrica.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp, symmetrica]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/symmetrica',
    install: true,
    include_directories: [],
    dependencies: deps,
  )
  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'libs/symmetrica',
      install: true,
      copy: true,
    )
  endif
endforeach

