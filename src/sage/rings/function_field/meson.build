py.install_sources(
  '__init__.py',
  'all.py',
  'constructor.py',
  'derivations.py',
  'derivations_polymod.py',
  'derivations_rational.py',
  'differential.py',
  'divisor.py',
  'element.pxd',
  'extensions.py',
  'function_field.py',
  'function_field_polymod.py',
  'function_field_rational.py',
  'ideal.py',
  'ideal_polymod.py',
  'ideal_rational.py',
  'jacobian_base.py',
  'jacobian_hess.py',
  'jacobian_khuri_makdisi.py',
  'maps.py',
  'order.py',
  'order_basis.py',
  'order_polymod.py',
  'order_rational.py',
  'place.py',
  'place_polymod.py',
  'place_rational.py',
  'valuation.py',
  'valuation_ring.py',
  subdir: 'sage/rings/function_field',
)

extension_data = {
  'element' : files('element.pyx'),
  'element_polymod' : files('element_polymod.pyx'),
  'element_rational' : files('element_rational.pyx'),
  'hermite_form_polynomial' : files('hermite_form_polynomial.pyx'),
  'khuri_makdisi' : files('khuri_makdisi.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/function_field',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp],
  )
endforeach

install_subdir(
  'drinfeld_modules',
  install_dir: sage_install_dir / 'rings/function_field',
)
