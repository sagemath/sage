py.install_sources(
  '__init__.py',
  'all.py',
  'gauss_legendre.pyx',
  'interactive_simplex_method.py',
  'knapsack.py',
  'linear_functions.pxd',
  'linear_functions.pyx',
  'linear_tensor.py',
  'linear_tensor_constraints.py',
  'linear_tensor_element.pxd',
  'linear_tensor_element.pyx',
  'mip.pxd',
  'mip.pyx',
  'optimize.py',
  'sdp.pxd',
  'sdp.pyx',
  subdir: 'sage/numerical',
)

extension_data = {
  'gauss_legendre' : files('gauss_legendre.pyx'),
  'linear_functions' : files('linear_functions.pyx'),
  'linear_tensor_element' : files('linear_tensor_element.pyx'),
  'mip' : files('mip.pyx'),
  'sdp' : files('sdp.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/numerical',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cypari2, gmp, mpfr],
  )
endforeach

subdir('backends')
