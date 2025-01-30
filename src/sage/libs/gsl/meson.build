py.install_sources(
  '__init__.py',
  'airy.pxd',
  'all.pxd',
  'array.pxd',
  'bessel.pxd',
  'blas.pxd',
  'blas_types.pxd',
  'block.pxd',
  'chebyshev.pxd',
  'clausen.pxd',
  'combination.pxd',
  'complex.pxd',
  'coulomb.pxd',
  'coupling.pxd',
  'dawson.pxd',
  'debye.pxd',
  'dilog.pxd',
  'eigen.pxd',
  'elementary.pxd',
  'ellint.pxd',
  'elljac.pxd',
  'erf.pxd',
  'errno.pxd',
  'exp.pxd',
  'expint.pxd',
  'fermi_dirac.pxd',
  'fft.pxd',
  'fit.pxd',
  'gamma.pxd',
  'gegenbauer.pxd',
  'histogram.pxd',
  'hyperg.pxd',
  'integration.pxd',
  'interp.pxd',
  'laguerre.pxd',
  'lambert.pxd',
  'legendre.pxd',
  'linalg.pxd',
  'log.pxd',
  'math.pxd',
  'matrix.pxd',
  'matrix_complex.pxd',
  'min.pxd',
  'monte.pxd',
  'ntuple.pxd',
  'odeiv.pxd',
  'permutation.pxd',
  'poly.pxd',
  'pow_int.pxd',
  'psi.pxd',
  'qrng.pxd',
  'random.pxd',
  'rng.pxd',
  'roots.pxd',
  'sort.pxd',
  'statistics.pxd',
  'sum.pxd',
  'synchrotron.pxd',
  'transport.pxd',
  'trig.pxd',
  'types.pxd',
  'vector.pxd',
  'vector_complex.pxd',
  'wavelet.pxd',
  'zeta.pxd',
  subdir: 'sage/libs/gsl',
)

extension_data = {'array' : files('array.pyx')}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp, gsl]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/gsl',
    install: true,
    include_directories: [inc_gsl],
    dependencies: deps,
  )

  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'libs/gsl',
      install: true,
      copy: true,
    )
  endif
endforeach

