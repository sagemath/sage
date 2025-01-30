iml = cc.find_library('iml')


py.install_sources(
  '__init__.py',
  'action.pxd',
  'all.py',
  'all__sagemath_meataxe.py',
  'args.pxd',
  'benchmark.py',
  'berlekamp_massey.py',
  'compute_J_ideal.py',
  'docs.py',
  'matrix.pxd',
  'matrix0.pxd',
  'matrix1.pxd',
  'matrix2.pxd',
  'matrix_cdv.pxd',
  'matrix_complex_ball_dense.pxd',
  'matrix_complex_double_dense.pxd',
  'matrix_cyclo_dense.pxd',
  'matrix_dense.pxd',
  'matrix_domain_dense.pxd',
  'matrix_domain_sparse.pxd',
  'matrix_double_dense.pxd',
  'matrix_double_sparse.pxd',
  'matrix_gap.pxd',
  'matrix_generic_dense.pxd',
  'matrix_generic_sparse.pxd',
  'matrix_gf2e_dense.pxd',
  'matrix_gfpn_dense.pxd',
  'matrix_integer_dense.pxd',
  'matrix_integer_dense_hnf.py',
  'matrix_integer_dense_saturation.py',
  'matrix_integer_sparse.pxd',
  'matrix_laurent_mpolynomial_dense.pxd',
  'matrix_misc.py',
  'matrix_mod2_dense.pxd',
  'matrix_modn_dense_double.pxd',
  'matrix_modn_dense_float.pxd',
  'matrix_modn_sparse.pxd',
  'matrix_mpolynomial_dense.pxd',
  'matrix_numpy_dense.pxd',
  'matrix_numpy_integer_dense.pxd',
  'matrix_polynomial_dense.pxd',
  'matrix_rational_dense.pxd',
  'matrix_rational_sparse.pxd',
  'matrix_real_double_dense.pxd',
  'matrix_space.py',
  'matrix_sparse.pxd',
  'matrix_symbolic_dense.pxd',
  'matrix_symbolic_sparse.pxd',
  'matrix_window.pxd',
  'operation_table.py',
  'special.py',
  'symplectic_basis.py',
  'template.pxd',
  'tests.py',
  subdir: 'sage/matrix',
)

extension_data = {
  'action' : files('action.pyx'),
  'args' : files('args.pyx'),
  'change_ring' : files('change_ring.pyx'),
  'constructor' : files('constructor.pyx'),
  'echelon_matrix' : files('echelon_matrix.pyx'),
  'matrix0' : files('matrix0.pyx'),
  'matrix1' : files('matrix1.pyx'),
  'matrix2' : files('matrix2.pyx'),
  'matrix_cdv' : files('matrix_cdv.pyx'),
  'matrix_complex_ball_dense' : files('matrix_complex_ball_dense.pyx'),
  'matrix_complex_double_dense' : files('matrix_complex_double_dense.pyx'),
  'matrix_dense' : files('matrix_dense.pyx'),
  'matrix_double_dense' : files('matrix_double_dense.pyx'),
  'matrix_double_sparse' : files('matrix_double_sparse.pyx'),
  'matrix_gap' : files('matrix_gap.pyx'),
  'matrix_generic_dense' : files('matrix_generic_dense.pyx'),
  'matrix_generic_sparse' : files('matrix_generic_sparse.pyx'),
  'matrix_gfpn_dense' : files('matrix_gfpn_dense.pyx'),
  'matrix_laurent_mpolynomial_dense' : files(
    'matrix_laurent_mpolynomial_dense.pyx',
  ),
  'matrix_numpy_dense' : files('matrix_numpy_dense.pyx'),
  'matrix_numpy_integer_dense' : files('matrix_numpy_integer_dense.pyx'),
  'matrix_polynomial_dense' : files('matrix_polynomial_dense.pyx'),
  'matrix_rational_sparse' : files('matrix_rational_sparse.pyx'),
  'matrix_real_double_dense' : files('matrix_real_double_dense.pyx'),
  'matrix_sparse' : files('matrix_sparse.pyx'),
  'matrix_symbolic_dense' : files('matrix_symbolic_dense.pyx'),
  'matrix_symbolic_sparse' : files('matrix_symbolic_sparse.pyx'),
  'matrix_window' : files('matrix_window.pyx'),
  'misc' : files('misc.pyx'),
  'misc_flint' : files('misc_flint.pyx'),
  'misc_mpfr' : files('misc_mpfr.pyx'),
  'strassen' : files('strassen.pyx'),
}

foreach name, pyx : extension_data
  dependencies = [py_dep, cysignals, gmp]
  if name == 'matrix_gfpn_dense'
    dependencies += [mtx]
  elif name == 'matrix_gap'
    dependencies += [gap]
  elif name == 'misc_mpfr'
    dependencies += [mpfr]
  elif name == 'matrix_complex_ball_dense'
    dependencies += [mpfi]
  elif name == 'misc_flint' or name == 'matrix_rational_sparse'
    dependencies += [flint]
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/matrix',
    install: true,
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_flint,
      inc_ntl,
      inc_numpy,
      inc_rings,
      inc_rings_finite,
    ],
    dependencies: dependencies,
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
      install_dir: sage_install_dir / 'matrix',
      install: true,
      copy: true,
    )
  endif
endforeach

extension_data_cpp = {
  'matrix_cyclo_dense': files('matrix_cyclo_dense.pyx'),
  'matrix_gf2e_dense': files('matrix_gf2e_dense.pyx'),
  'matrix_integer_dense': files('matrix_integer_dense.pyx'),
  'matrix_integer_sparse': files('matrix_integer_sparse.pyx'),
  'matrix_mod2_dense': files('matrix_mod2_dense.pyx'),
  'matrix_modn_dense_double': files('matrix_modn_dense_double.pyx'),
  'matrix_modn_dense_float': files('matrix_modn_dense_float.pyx'),
  'matrix_modn_sparse': files('matrix_modn_sparse.pyx'),
  'matrix_mpolynomial_dense': files('matrix_mpolynomial_dense.pyx'),
  'matrix_rational_dense': files('matrix_rational_dense.pyx'),
}

foreach name, pyx : extension_data_cpp
  deps = [
    py_dep,
    blas,
    cypari2,
    cysignals,
    fflas,
    flint,
    gd,
    givaro,
    gmp,
    gmpxx,
    iml,
    linbox,
    m,
    m4ri,
    m4rie,
    mpfi,
    mpfr,
    ntl,
    pari,
    png,
    singular,
    zlib,
  ]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/matrix',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_flint,
      inc_ntl,
      inc_numpy,
      inc_rings,
      inc_rings_finite,
    ],
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
      install_dir: sage_install_dir / 'matrix',
      install: true,
      copy: true,
    )
  endif
endforeach

