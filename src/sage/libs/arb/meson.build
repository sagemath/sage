py.install_sources(
  '__init__.py',
  'arb_wrap.h',
  'arith.pyx',
  subdir: 'sage/libs/arb',
)

extension_data = {'arith': files('arith.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, flint, gmp, mpfr]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/arb',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings],
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
      install_dir: sage_install_dir / 'libs/arb',
      install: true,
      copy: true,
    )
  endif
endforeach

