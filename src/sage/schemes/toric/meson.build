py.install_sources(
  'all.py',
  'chow_group.py',
  'divisor.py',
  'fano_variety.py',
  'homset.py',
  'ideal.py',
  'library.py',
  'morphism.py',
  'points.py',
  'toric_subscheme.py',
  'variety.py',
  'weierstrass.py',
  'weierstrass_covering.py',
  'weierstrass_higher.py',
  subdir: 'sage/schemes/toric',
)

extension_data = {'divisor_class' : files('divisor_class.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/schemes/toric',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp],
  )
endforeach

install_subdir('sheaf', install_dir: sage_install_dir / 'schemes/toric')
