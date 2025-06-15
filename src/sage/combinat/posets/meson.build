py.install_sources(
  '__init__.py',
  'all.py',
  'bubble_shuffle.py',
  'cartesian_product.py',
  'd_complete.py',
  'elements.py',
  'forest.py',
  'hasse_diagram.py',
  'hochschild_lattice.py',
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
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/posets',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

