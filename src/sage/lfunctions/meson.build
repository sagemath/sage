py.install_sources(
  '__init__.py',
  'all.py',
  'dokchitser.py',
  'lcalc.py',
  'pari.py',
  'sympow.py',
  subdir: 'sage/lfunctions',
)

extension_data = {'zero_sums' : files('zero_sums.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, flint, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/lfunctions',
    install: true,
    include_directories: [inc_flint],
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
      install_dir: sage_install_dir / 'lfunctions',
      install: true,
      copy: true,
    )
  endif
endforeach

