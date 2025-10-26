# Cannot be found via pkg-config
glpk = cc.find_library('glpk', required: not is_windows, disabler: true)
highs = dependency('highs', required: false, disabler: true)
if not highs.found()
  highs = cc.find_library('highs', required: not is_windows, disabler: true)
endif
# HiGHS headers use relative includes like #include "lp_data/HighsCallbackStruct.h"
# so we need to add the highs subdirectory to the include path.
# Check multiple locations in order: SAGE_LOCAL, Python prefix (conda), system paths
highs_compile_args = []

# 1. Check SAGE_LOCAL
sage_local = get_option('SAGE_LOCAL')
if sage_local == ''
  sage_local = meson.project_source_root() / '..' / 'local'
endif
highs_local_include = sage_local / 'include' / 'highs'
if fs.is_dir(highs_local_include)
  highs_compile_args += ['-I' + highs_local_include]
endif

# 2. Check Python prefix (for conda environments) - this is the most common case for CI
python_prefix = run_command(
  py,
  ['-c', 'import sys; print(sys.prefix)'],
  check: false
).stdout().strip()

if python_prefix != ''
  highs_prefix_include = python_prefix / 'include' / 'highs'
  if fs.is_dir(highs_prefix_include)
    highs_compile_args += ['-I' + highs_prefix_include]
  endif
endif

# 3. Check common system paths
foreach include_base : ['/usr/include', '/usr/local/include', get_option('prefix') / 'include']
  highs_system_include = include_base / 'highs'
  if fs.is_dir(highs_system_include)
    highs_compile_args += ['-I' + highs_system_include]
    break
  endif
endforeach
py.install_sources(
  '__init__.py',
  'all.py',
  'cvxopt_backend.pyx',
  'cvxopt_backend_test.py',
  'cvxopt_sdp_backend.pyx',
  'cvxpy_backend.pxd',
  'cvxpy_backend.pyx',
  'cvxpy_backend_test.py',
  'generic_backend.pxd',
  'generic_backend.pyx',
  'generic_backend_test.py',
  'generic_sdp_backend.pxd',
  'generic_sdp_backend.pyx',
  'glpk_backend.pxd',
  'glpk_backend.pyx',
  'glpk_backend_test.py',
  'glpk_exact_backend.pxd',
  'glpk_exact_backend.pyx',
  'glpk_exact_backend_test.py',
  'glpk_graph_backend.pxd',
  'glpk_graph_backend.pyx',
  'highs_backend.pxd',
  'highs_backend.pyx',
  'interactivelp_backend.pxd',
  'interactivelp_backend.pyx',
  'interactivelp_backend_test.py',
  'logging_backend.py',
  'matrix_sdp_backend.pxd',
  'matrix_sdp_backend.pyx',
  'ppl_backend.pyx',
  'ppl_backend_test.py',
  'scip_backend.pxd',
  'scip_backend.pyx',
  'scip_backend_test.py',
  subdir: 'sage/numerical/backends',
)

extension_data = {
  'cvxopt_backend' : files('cvxopt_backend.pyx'),
  'cvxopt_sdp_backend' : files('cvxopt_sdp_backend.pyx'),
  'cvxpy_backend' : files('cvxpy_backend.pyx'),
  'generic_backend' : files('generic_backend.pyx'),
  'generic_sdp_backend' : files('generic_sdp_backend.pyx'),
  'glpk_backend' : files('glpk_backend.pyx'),
  'glpk_exact_backend' : files('glpk_exact_backend.pyx'),
  'glpk_graph_backend' : files('glpk_graph_backend.pyx'),
  'highs_backend' : files('highs_backend.pyx'),
  'interactivelp_backend' : files('interactivelp_backend.pyx'),
  'matrix_sdp_backend' : files('matrix_sdp_backend.pyx'),
  'ppl_backend' : files('ppl_backend.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/numerical/backends',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    c_args: highs_compile_args,
    dependencies: [py_dep, cysignals, glpk, highs, gmp],
  )
endforeach

extension_data_cpp = {'scip_backend': files('scip_backend.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/numerical/backends',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach
