lcalc = dependency('lcalc', version: '>= 2.0.0')

py.install_sources(
  '__init__.py',
  'lcalc_Lfunction.pxd',
  'lcalc_sage.h',
  subdir: 'sage/libs/lcalc',
)

extension_data_cpp = {'lcalc_Lfunction': files('lcalc_Lfunction.pyx')}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cypari2, cysignals, gmp, lcalc, m, mpfr, ntl]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/lcalc',
    install: true,
    override_options: ['cython_language=cpp'],
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
      install_dir: sage_install_dir / 'libs/lcalc',
      install: true,
      copy: true,
    )
  endif
endforeach

