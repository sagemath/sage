py.install_sources(
  'abelian.py',
  'affine_lie_algebra.py',
  'all.py',
  'bch.py',
  'bgg_dual_module.py',
  'bgg_resolution.py',
  'center_uea.py',
  'classical_lie_algebra.py',
  'examples.py',
  'free_lie_algebra.py',
  'heisenberg.py',
  'lie_algebra.py',
  'lie_algebra_element.pxd',
  'morphism.py',
  'nilpotent_lie_algebra.py',
  'onsager.py',
  'poincare_birkhoff_witt.py',
  'quotient.py',
  'rank_two_heisenberg_virasoro.py',
  'representation.py',
  'structure_coefficients.py',
  'subalgebra.py',
  'symplectic_derivation.py',
  'verma_module.py',
  'virasoro.py',
  subdir: 'sage/algebras/lie_algebras',
)

extension_data = {'lie_algebra_element' : files('lie_algebra_element.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/lie_algebras',
    install: true,
    include_directories: [inc_cpython, inc_data_structures],
    dependencies: [py_dep, gmp],
  )
endforeach

