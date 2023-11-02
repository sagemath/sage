# Cannot be found via pkg-config
meataxe = declare_dependency(
  dependencies: [
    cc.find_library('meataxe', required: false),
  ]
)
iml = declare_dependency(
  dependencies: [
    cc.find_library('iml'),
  ]
)


py.install_sources(
    'all.py',
    'all__sagemath_meataxe.py',
    'benchmark.py',
    'berlekamp_massey.py',
    'compute_J_ideal.py',
    'docs.py',
    'matrix_integer_dense_hnf.py',
    'matrix_integer_dense_saturation.py',
    'matrix_misc.py',
    'matrix_space.py',
    'operation_table.py',
    'special.py',
    'symplectic_basis.py',
    'tests.py',
    subdir: 'sage/matrix',
)


extension_data = {
    'action': files('action.pyx'),
    'args': files('args.pyx'),
    'change_ring': files('change_ring.pyx'),
    'constructor': files('constructor.pyx'),
    'echelon_matrix': files('echelon_matrix.pyx'),
    'matrix0': files('matrix0.pyx'),
    'matrix1': files('matrix1.pyx'),
    'matrix2': files('matrix2.pyx'),
    'matrix_cdv': files('matrix_cdv.pyx'),
    'matrix_complex_ball_dense': files('matrix_complex_ball_dense.pyx'),
    'matrix_complex_double_dense': files('matrix_complex_double_dense.pyx'),
    'matrix_dense': files('matrix_dense.pyx'),
    'matrix_double_dense': files('matrix_double_dense.pyx'),
    'matrix_double_sparse': files('matrix_double_sparse.pyx'),
    'matrix_gap': files('matrix_gap.pyx'),
    'matrix_generic_dense': files('matrix_generic_dense.pyx'),
    'matrix_generic_sparse': files('matrix_generic_sparse.pyx'),
    'matrix_gf2e_dense': files('matrix_gf2e_dense.pyx'),
    #'matrix_gfpn_dense': files('matrix_gfpn_dense.pyx'), - requires meataxe
    'matrix_integer_dense': files('matrix_integer_dense.pyx'),
    'matrix_mod2_dense': files('matrix_mod2_dense.pyx'),
    #'matrix_mpolynomial_dense': files('matrix_mpolynomial_dense.pyx'), - compilation error
    'matrix_numpy_dense': files('matrix_numpy_dense.pyx'),
    'matrix_numpy_integer_dense': files('matrix_numpy_integer_dense.pyx'),
    'matrix_polynomial_dense': files('matrix_polynomial_dense.pyx'),
    'matrix_rational_sparse': files('matrix_rational_sparse.pyx'),
    'matrix_real_double_dense': files('matrix_real_double_dense.pyx'),
    'matrix_sparse': files('matrix_sparse.pyx'),
    'matrix_symbolic_dense': files('matrix_symbolic_dense.pyx'),
    'matrix_symbolic_sparse': files('matrix_symbolic_sparse.pyx'),
    'matrix_window': files('matrix_window.pyx'),
    'misc': files('misc.pyx'),
    'misc_flint': files('misc_flint.pyx'),
    'misc_mpfr': files('misc_mpfr.pyx'),
    'strassen': files('strassen.pyx'),
}
extension_data_cpp = {
    'matrix_cyclo_dense': files('matrix_cyclo_dense.pyx'),
    'matrix_modn_dense_double': files('matrix_modn_dense_double.pyx'),
    'matrix_modn_dense_float': files('matrix_modn_dense_float.pyx'),
    'matrix_integer_sparse': files('matrix_integer_sparse.pyx'),
    'matrix_modn_sparse': files('matrix_modn_sparse.pyx'),
    'matrix_rational_dense': files('matrix_rational_dense.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/matrix',
        install: true,
        include_directories: [inc_cpython, inc_numpy, inc_flint, inc_ext, inc_rings, inc_rings_finite, inc_arb, inc_ntl, inc_iml],
        dependencies: [py_dep, cysignals, cypari2, gmp, mpfr, pari, flint, iml],
    )
endforeach
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/matrix',
        install: true,
        override_options: ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_numpy, inc_flint, inc_ext, inc_rings_finite, inc_rings, inc_arb, inc_ntl],
        dependencies: [py_dep, cysignals, cypari2, gmp, mpfr],
    )
endforeach
