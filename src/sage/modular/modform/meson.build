py.install_sources(
  'all.py',
  'ambient.py',
  'ambient_R.py',
  'ambient_eps.py',
  'ambient_g0.py',
  'ambient_g1.py',
  'constructor.py',
  'cuspidal_submodule.py',
  'defaults.py',
  'eis_series.py',
  'eisenstein_submodule.py',
  'element.py',
  'find_generators.py',
  'half_integral.py',
  'hecke_operator_on_qexp.py',
  'j_invariant.py',
  'l_series_gross_zagier.py',
  'notes.py',
  'numerical.py',
  'periods.py',
  'ring.py',
  'space.py',
  'submodule.py',
  'tests.py',
  'theta.py',
  'vm_basis.py',
  'weight1.py',
  subdir: 'sage/modular/modform',
)

extension_data = {
  'eis_series_cython' : files('eis_series_cython.pyx'),
  'l_series_gross_zagier_coeffs' : files('l_series_gross_zagier_coeffs.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular/modform',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

