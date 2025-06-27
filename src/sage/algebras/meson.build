py.install_sources(
  '__init__.py',
  'affine_nil_temperley_lieb.py',
  'all.py',
  'askey_wilson.py',
  'associated_graded.py',
  'catalog.py',
  'cellular_basis.py',
  'clifford_algebra.py',
  'clifford_algebra_element.pxd',
  'cluster_algebra.py',
  'commutative_dga.py',
  'down_up_algebra.py',
  'exterior_algebra_groebner.pxd',
  'finite_gca.py',
  'free_algebra.py',
  'free_algebra_element.py',
  'free_algebra_quotient.py',
  'free_algebra_quotient_element.py',
  'free_zinbiel_algebra.py',
  'group_algebra.py',
  'hall_algebra.py',
  'iwahori_hecke_algebra.py',
  'jordan_algebra.py',
  'nil_coxeter_algebra.py',
  'octonion_algebra.pxd',
  'orlik_solomon.py',
  'orlik_terao.py',
  'q_commuting_polynomials.py',
  'q_system.py',
  'quantum_clifford.py',
  'quantum_matrix_coordinate_algebra.py',
  'quantum_oscillator.py',
  'quaternion_algebra.py',
  'quaternion_algebra_element.py',
  'rational_cherednik_algebra.py',
  'schur_algebra.py',
  'shuffle_algebra.py',
  'splitting_algebra.py',
  'tensor_algebra.py',
  'weyl_algebra.py',
  'yangian.py',
  'yokonuma_hecke_algebra.py',
  subdir: 'sage/algebras',
)

extension_data = {
  'clifford_algebra_element' : files('clifford_algebra_element.pyx'),
  'exterior_algebra_groebner' : files('exterior_algebra_groebner.pyx'),
  'octonion_algebra' : files('octonion_algebra.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

subdir('finite_dimensional_algebras')
subdir('fusion_rings')
install_subdir('hecke_algebras', install_dir: sage_install_dir / 'algebras')
subdir('letterplace')
subdir('lie_algebras')
install_subdir(
  'lie_conformal_algebras',
  install_dir: sage_install_dir / 'algebras',
)
install_subdir('quantum_groups', install_dir: sage_install_dir / 'algebras')
subdir('quatalg')
install_subdir('steenrod', install_dir: sage_install_dir / 'algebras')
