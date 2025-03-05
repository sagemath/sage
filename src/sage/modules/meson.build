py.install_sources(
  '__init__.py',
  'all.py',
  'complex_double_vector.py',
  'diamond_cutting.py',
  'filtered_vector_space.py',
  'finite_submodule_iter.pxd',
  'free_module.py',
  'free_module_element.pxd',
  'free_module_homspace.py',
  'free_module_integer.py',
  'free_module_morphism.py',
  'free_module_pseudohomspace.py',
  'free_module_pseudomorphism.py',
  'free_quadratic_module.py',
  'free_quadratic_module_integer_symmetric.py',
  'matrix_morphism.py',
  'misc.py',
  'module.pxd',
  'module_functors.py',
  'multi_filtered_vector_space.py',
  'numpy_util.pxd',
  'quotient_module.py',
  'real_double_vector.py',
  'submodule.py',
  'tensor_operations.py',
  'torsion_quadratic_module.py',
  'tutorial_free_modules.py',
  'vector_callable_symbolic_dense.py',
  'vector_complex_double_dense.pxd',
  'vector_double_dense.pxd',
  'vector_integer_dense.pxd',
  'vector_integer_sparse.pxd',
  'vector_mod2_dense.pxd',
  'vector_modn_dense.pxd',
  'vector_modn_sparse.pxd',
  'vector_numpy_dense.pxd',
  'vector_numpy_integer_dense.pxd',
  'vector_rational_dense.pxd',
  'vector_rational_sparse.pxd',
  'vector_real_double_dense.pxd',
  'vector_space_homspace.py',
  'vector_space_morphism.py',
  'vector_symbolic_dense.py',
  'vector_symbolic_sparse.py',
  subdir: 'sage/modules',
)

extension_data = {
  'finite_submodule_iter' : files('finite_submodule_iter.pyx'),
  'free_module_element' : files('free_module_element.pyx'),
  'module' : files('module.pyx'),
  'vector_complex_double_dense' : files('vector_complex_double_dense.pyx'),
  'vector_double_dense' : files('vector_double_dense.pyx'),
  'vector_integer_dense' : files('vector_integer_dense.pyx'),
  'vector_integer_sparse' : files('vector_integer_sparse.pyx'),
  'vector_modn_dense' : files('vector_modn_dense.pyx'),
  'vector_modn_sparse' : files('vector_modn_sparse.pyx'),
  'vector_numpy_dense' : files('vector_numpy_dense.pyx'),
  'vector_numpy_integer_dense' : files('vector_numpy_integer_dense.pyx'),
  'vector_rational_dense' : files('vector_rational_dense.pyx'),
  'vector_rational_sparse' : files('vector_rational_sparse.pyx'),
  'vector_real_double_dense' : files('vector_real_double_dense.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modules',
    install: true,
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_numpy,
      inc_rings,
      inc_rings_finite,
    ],
    dependencies: [py_dep, cysignals, gd, gmp, m4ri, png],
  )
endforeach

extension_data_cpp = {
  'numpy_util' : files('numpy_util.pyx'),
  'vector_mod2_dense': files('vector_mod2_dense.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modules',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_numpy,
      inc_rings,
      inc_rings_finite,
    ],
    dependencies: [py_dep, cysignals, gd, gmp, m4ri, png],
  )
endforeach

install_subdir('fg_pid', install_dir: sage_install_dir / 'modules')
install_subdir('fp_graded', install_dir: sage_install_dir / 'modules')
subdir('with_basis')
