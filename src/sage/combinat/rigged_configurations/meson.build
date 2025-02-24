py.install_sources(
  'all.py',
  'bij_abstract_class.py',
  'bij_infinity.py',
  'bij_type_A.py',
  'bij_type_A2_dual.py',
  'bij_type_A2_even.py',
  'bij_type_A2_odd.py',
  'bij_type_B.py',
  'bij_type_C.py',
  'bij_type_D.py',
  'bij_type_D_tri.py',
  'bij_type_D_twisted.py',
  'bij_type_E67.py',
  'bijection.py',
  'kleber_tree.py',
  'kr_tableaux.py',
  'rc_crystal.py',
  'rc_infinity.py',
  'rigged_configuration_element.py',
  'rigged_configurations.py',
  'rigged_partition.pxd',
  'tensor_product_kr_tableaux.py',
  'tensor_product_kr_tableaux_element.py',
  subdir: 'sage/combinat/rigged_configurations',
)

extension_data = {'rigged_partition' : files('rigged_partition.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/rigged_configurations',
    install: true,
    include_directories: [],
    dependencies: [py_dep, gmp],
  )
endforeach

