py.install_sources(
  '__init__.py',
  'affine_ds.py',
  'all.py',
  'berkovich_ds.py',
  'dynamical_semigroup.py',
  'endPN_automorphism_group.py',
  'endPN_minimal_model.py',
  'generic_ds.py',
  'product_projective_ds.py',
  'projective_ds.py',
  'wehlerK3.py',
  subdir: 'sage/dynamics/arithmetic_dynamics',
)

extension_data = {'projective_ds_helper' : files('projective_ds_helper.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/dynamics/arithmetic_dynamics',
    install: true,
    include_directories: [],
    dependencies: [py_dep, gmp],
  )
endforeach

