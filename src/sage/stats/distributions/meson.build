py.install_sources(
  '__init__.py',
  'all.py',
  'catalog.py',
  'dgs.pxd',
  'discrete_gaussian_integer.pxd',
  'discrete_gaussian_lattice.py',
  'discrete_gaussian_polynomial.py',
  subdir: 'sage/stats/distributions',
)

extension_data = {
  'discrete_gaussian_integer' : files(
    'dgs_bern.c',
    'dgs_gauss_dp.c',
    'dgs_gauss_mp.c',
    'discrete_gaussian_integer.pyx',
  ),
}

foreach name, pyx : extension_data
  deps = [py_dep, cypari2, cysignals, gmp, mpfr]
  py.extension_module(
    name,
    sources: pyx,
    c_args: ['-D_XOPEN_SOURCE=600'],
    subdir: 'sage/stats/distributions',
    install: true,
    include_directories: [inc_cpython, inc_rings, inc_src],
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
      install_dir: sage_install_dir / 'stats/distributions',
      install: true,
      copy: true,
    )
  endif
endforeach

