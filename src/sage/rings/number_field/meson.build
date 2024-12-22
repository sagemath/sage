py.install_sources(
  'S_unit_solver.py',
  '__init__.py',
  'all.py',
  'bdd_height.py',
  'class_group.py',
  'galois_group.py',
  'homset.py',
  'maps.py',
  'morphism.py',
  'number_field.py',
  'number_field_base.pxd',
  'number_field_element.pxd',
  'number_field_element_base.pxd',
  'number_field_element_quadratic.pxd',
  'number_field_ideal.py',
  'number_field_ideal_rel.py',
  'number_field_rel.py',
  'order.py',
  'order_ideal.py',
  'selmer_group.py',
  'small_primes_of_degree_one.py',
  'splitting_field.py',
  'structure.py',
  'totallyreal_data.pxd',
  'totallyreal_phc.py',
  'totallyreal_rel.py',
  'unit_group.py',
  subdir: 'sage/rings/number_field',
)

extension_data = {
  'number_field_base' : files('number_field_base.pyx'),
  'number_field_element_base' : files('number_field_element_base.pyx'),
  'number_field_morphisms' : files('number_field_morphisms.pyx'),
  'totallyreal' : files('totallyreal.pyx'),
  'totallyreal_data' : files('totallyreal_data.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name == 'totallyreal'
    deps += [cypari2]
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/number_field',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_ntl, inc_rings],
    dependencies: deps,
  )
endforeach

extension_data_cpp = {
  'number_field_element': files('number_field_element.pyx'),
  'number_field_element_quadratic': files('number_field_element_quadratic.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/number_field',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_ntl, inc_rings],
    dependencies: [py_dep, cypari2, cysignals, flint, gmp, mpfi, mpfr, ntl],
  )
endforeach

