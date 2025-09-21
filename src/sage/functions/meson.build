py.install_sources(
  '__init__.py',
  'airy.py',
  'all.py',
  'bessel.py',
  'error.py',
  'exp_integral.py',
  'gamma.py',
  'generalized.py',
  'hyperbolic.py',
  'hypergeometric.py',
  'jacobi.py',
  'log.py',
  'min_max.py',
  'orthogonal_polys.py',
  'other.py',
  'piecewise.py',
  'special.py',
  'spike_function.py',
  'transcendental.py',
  'trig.py',
  'wigner.py',
  subdir: 'sage/functions',
)

extension_data = {'prime_pi' : files('prime_pi.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/functions',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp],
  )
endforeach

