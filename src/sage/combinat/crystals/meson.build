py.install_sources(
  'affine.py',
  'affine_factorization.py',
  'affinization.py',
  'alcove_path.py',
  'all.py',
  'bkk_crystals.py',
  'catalog.py',
  'catalog_elementary_crystals.py',
  'catalog_infinity_crystals.py',
  'catalog_kirillov_reshetikhin.py',
  'crystals.py',
  'direct_sum.py',
  'elementary_crystals.py',
  'fast_crystals.py',
  'fully_commutative_stable_grothendieck.py',
  'generalized_young_walls.py',
  'highest_weight_crystals.py',
  'induced_structure.py',
  'infinity_crystals.py',
  'kac_modules.py',
  'kirillov_reshetikhin.py',
  'kyoto_path_model.py',
  'letters.pxd',
  'littelmann_path.py',
  'monomial_crystals.py',
  'multisegments.py',
  'mv_polytopes.py',
  'pbw_crystal.py',
  'pbw_datum.pxd',
  'polyhedral_realization.py',
  'spins.pxd',
  'star_crystal.py',
  'subcrystal.py',
  'tensor_product.py',
  'tensor_product_element.pxd',
  'virtual_crystal.py',
  subdir: 'sage/combinat/crystals',
)

extension_data = {
  'letters' : files('letters.pyx'),
  'pbw_datum' : files('pbw_datum.pyx'),
  'spins' : files('spins.pyx'),
  'tensor_product_element' : files('tensor_product_element.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/crystals',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

