coxeter3 = cc.find_library('coxeter3', required: false, disabler: true)
py.install_sources(
  '__init__.py',
  'all__sagemath_coxeter3.py',
  'coxeter.pxd',
  'coxeter_group.py',
  'decl.pxd',
  subdir: 'sage/libs/coxeter3',
)

extension_data_cpp = {'coxeter': files('coxeter.pyx')}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cysignals, coxeter3]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/coxeter3',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython],
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
      install_dir: sage_install_dir / 'libs/coxeter3',
      install: true,
      copy: true,
    )
  endif
endforeach
