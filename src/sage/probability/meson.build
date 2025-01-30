py.install_sources(
  '__init__.py',
  'all.py',
  'random_variable.py',
  subdir: 'sage/probability',
)

extension_data = {
  'probability_distribution' : files('probability_distribution.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp, gsl]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/probability',
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
      install_dir: sage_install_dir / 'probability',
      install: true,
      copy: true,
    )
  endif
endforeach

