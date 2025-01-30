py.install_sources(
  '__init__.py',
  'all.py',
  'mpfi.pxd',
  subdir: 'sage/rings/convert',
)

extension_data = {'mpfi' : files('mpfi.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, gmp, gsl, mpfi, mpfr]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/convert',
    install: true,
    include_directories: [inc_cpython, inc_rings],
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
      install_dir: sage_install_dir / 'rings/convert',
      install: true,
      copy: true,
    )
  endif
endforeach

