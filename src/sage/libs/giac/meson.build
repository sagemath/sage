giac = cc.find_library('giac', required: false, disabler: true)

py.install_sources('__init__.py', 'giac.pxd', 'misc.h', subdir: 'sage/libs/giac')

extension_data_cpp = {'giac': files('giac.pyx')}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cysignals, giac, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/giac',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_rings],
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
      install_dir: sage_install_dir / 'libs/giac',
      install: true,
      copy: true,
    )
  endif
endforeach
