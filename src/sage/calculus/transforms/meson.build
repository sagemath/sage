py.install_sources(
  '__init__.py',
  'all.py',
  'dft.py',
  'dwt.pxd',
  'fft.pxd',
  subdir: 'sage/calculus/transforms',
)

extension_data = {'dwt' : files('dwt.pyx'), 'fft' : files('fft.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp, gsl]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/calculus/transforms',
    install: true,
    include_directories: [inc_gsl],
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
      install_dir: sage_install_dir / 'calculus/transforms',
      install: true,
      copy: true,
    )
  endif
endforeach

