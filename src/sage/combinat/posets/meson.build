py.install_sources(
  '__init__.py',
  'all.py',
  'cartesian_product.py',
  'd_complete.py',
  'elements.py',
  'forest.py',
  'hasse_diagram.py',
  'incidence_algebras.py',
  'lattices.py',
  'linear_extensions.py',
  'mobile.py',
  'moebius_algebra.py',
  'poset_examples.py',
  'posets.py',
  subdir: 'sage/combinat/posets',
)

extension_data = {
  'hasse_cython' : files('hasse_cython.pyx'),
  'hasse_cython_flint' : files('hasse_cython_flint.pyx'),
  'linear_extension_iterator' : files('linear_extension_iterator.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, flint, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/posets',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings],
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
      install_dir: sage_install_dir / 'combinat/posets',
      install: true,
      copy: true,
    )
  endif
endforeach

