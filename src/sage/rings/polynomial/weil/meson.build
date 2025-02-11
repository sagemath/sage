py.install_sources(
  '__init__.py',
  'all.py',
  'power_sums.h',
  subdir: 'sage/rings/polynomial/weil',
)

extension_data = {'weil_polynomials' : files('weil_polynomials.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, flint, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/polynomial/weil',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings, inc_src],
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
      install_dir: sage_install_dir / 'rings/polynomial/weil',
      install: true,
      copy: true,
    )
  endif
endforeach

