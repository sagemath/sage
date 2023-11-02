# Setup dependencies that are needed by many modules

inc_numpy = run_command(py,
    [
    '-c',
    '''
import numpy
print(numpy.get_include())
    '''.strip()
    ],
    check: true
  ).stdout().strip()
numpy = declare_dependency(
    include_directories: inc_numpy,
)

inc_cysignals = run_command(py,
    [
    '-c',
    '''
import cysignals
print(cysignals.__file__.replace('__init__.py', ''))
    '''.strip()
    ],
    check: true
  ).stdout().strip()
cysignals = declare_dependency(
    include_directories: inc_cysignals,
)

inc_gmpy2 = run_command(py,
    [
    '-c',
    '''
import gmpy2
print(gmpy2.__file__.replace('__init__.py', ''))
    '''.strip()
    ],
    check: true
  ).stdout().strip()
gmpy2 = declare_dependency(
    include_directories: inc_gmpy2,
)
gmp = dependency('gmp', required: true)

inc_cypari2 = run_command(py,
    [
    '-c',
    '''
import cypari2
print(cypari2.__file__.replace('__init__.py', ''))
    '''.strip()
    ],
    check: true
  ).stdout().strip()
cypari2 = declare_dependency(
    include_directories: inc_cypari2,
)
pari = declare_dependency(
  dependencies: [
    cc.find_library('pari'),
  ]
)

mpfr = declare_dependency(
  dependencies: [
    cc.find_library('mpfr'),
  ]
)

# Cannot be found by pkg-config
flint = declare_dependency(
  dependencies: [
    cc.find_library('flint'),
  ]
)

# It's strange but cython cannot find its own include files
# so we find them ourselves, and add them to the include path
inc_cython = run_command(py,
    [
    '-c',
    '''
import Cython
print(Cython.__file__.replace('__init__.py', ''))
    '''.strip()
    ],
    check: true
   ).stdout().strip()
add_project_arguments('-I', inc_cython + 'Includes', language : 'cython')

# Meson currently ignores include_directories for Cython modules, so we
# have to add them manually.
# https://github.com/mesonbuild/meson/issues/9562
add_project_arguments('-I', meson.current_source_dir(), language : 'cython')
#add_project_arguments('-I', py.get_path('include'), language : 'cython')

# Add global compiler flags
add_project_arguments('-X auto_pickle=False', language : 'cython')
add_project_arguments('-X autotestdict=False', language : 'cython')
add_project_arguments('-X binding=False', language : 'cython')
add_project_arguments('-X c_api_binop_methods=True', language : 'cython')
add_project_arguments('-X cdivision=True', language : 'cython')
add_project_arguments('-X cpow=True', language : 'cython')
add_project_arguments('-X embedsignature=True', language : 'cython')
add_project_arguments('-X fast_getattr=True', language : 'cython')
#add_project_arguments('-X language_level="3"', language : 'cython')
add_project_arguments('-X legacy_implicit_noexcept=True', language : 'cython')
add_project_arguments('-X preliminary_late_includes_cy28=True', language : 'cython')

inc_cpython = include_directories('sage/cpython')
inc_rings = include_directories('sage/rings')
inc_rings_finite = include_directories('sage/rings/finite_rings')
inc_flint = include_directories('sage/libs/flint')
inc_gsl = include_directories('sage/libs/gsl')
inc_ntl = include_directories('sage/libs/ntl')
inc_arb = include_directories('sage/libs/arb')
inc_iml = include_directories('sage/libs')
inc_data_structures = include_directories('sage/data_structures')
inc_ext = include_directories('sage/ext')
inc_src = include_directories('.')

# Submodules
subdir('sage')
