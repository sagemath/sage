py.install_sources(
  '__init__.py',
  'algebras_with_basis.py',
  'all.py',
  'commutative_additive_monoids.py',
  'commutative_additive_semigroups.py',
  'coxeter_groups.py',
  'crystals.py',
  'cw_complexes.py',
  'facade_sets.py',
  'filtered_algebras_with_basis.py',
  'filtered_modules_with_basis.py',
  'finite_coxeter_groups.py',
  'finite_dimensional_algebras_with_basis.py',
  'finite_dimensional_lie_algebras_with_basis.py',
  'finite_enumerated_sets.py',
  'finite_monoids.py',
  'finite_semigroups.py',
  'finite_weyl_groups.py',
  'graded_connected_hopf_algebras_with_basis.py',
  'graded_modules_with_basis.py',
  'graphs.py',
  'hopf_algebras_with_basis.py',
  'infinite_enumerated_sets.py',
  'lie_algebras.py',
  'lie_algebras_with_basis.py',
  'magmas.py',
  'manifolds.py',
  'monoids.py',
  'posets.py',
  'semigroups.py',
  'sets_cat.py',
  'sets_with_grading.py',
  'with_realizations.py',
  subdir: 'sage/categories/examples',
)

extension_data = {'semigroups_cython' : files('semigroups_cython.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/categories/examples',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, gmp],
  )
endforeach

