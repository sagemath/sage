# Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cython = meson.get_compiler('cython')

# Setup dependencies that are needed by many modules
numpy = dependency('numpy')

inc_cysignals = run_command(
  py,
  [
    '-c',
    '''
import cysignals
print(cysignals.__file__.replace('__init__.py', ''))
    '''.strip(),
  ],
  check: true,
).stdout().strip()
cysignals = declare_dependency(include_directories: inc_cysignals)

inc_gmpy2 = run_command(
  py,
  [
    '-c',
    '''
import gmpy2
print(gmpy2.__file__.replace('__init__.py', ''))
    '''.strip(),
  ],
  check: true,
).stdout().strip()
gmpy2 = declare_dependency(include_directories: inc_gmpy2)
gmp = dependency('gmp', required: true)

inc_cypari2 = run_command(
  py,
  [
    '-c',
    '''
import cypari2
print(cypari2.__file__.replace('__init__.py', ''))
    '''.strip(),
  ],
  check: true,
).stdout().strip()
cypari2 = declare_dependency(include_directories: inc_cypari2)
pari = cc.find_library('pari')

mpfr = cc.find_library('mpfr')

# Once https://github.com/mesonbuild/meson/pull/12616 is released, we can simply do the following:
# flint = dependency('flint', version : '>=3.0.0')
flint = dependency(
  'flint',
  version: '>=3.0.0',
  required: false,
  allow_fallback: false,
)
if not flint.found()
  # Currently the flint pkg-config file is broken, so we manually use find_library
  # Can be removed once we require a version of flint that ships https://github.com/flintlib/flint/pull/1647
  flint = cc.find_library('flint')
endif

blas_order = []
if host_machine.system() == 'darwin'
  blas_order += 'accelerate'
endif
if host_machine.cpu_family() == 'x86_64'
  blas_order += 'mkl'
endif
# pkg-config uses a lower-case name while CMake uses a capitalized name, so try
# that too to make the fallback detection with CMake work
blas_order += ['cblas', 'openblas', 'OpenBLAS', 'flexiblas', 'blis', 'blas']
blas = dependency(blas_order)
gsl = dependency(
  'gsl',
  fallback: ['gsl', 'gsl_dep'],
  version: '>=2.5',
  required: true,
)
gd = cc.find_library('gd')
m = cc.find_library('m')
m4ri = cc.find_library('m4ri')
m4rie = cc.find_library('m4rie')
mtx = cc.find_library('mtx', required: false)
png = cc.find_library('png')
zlib = cc.find_library('z')
intl = cc.find_library('intl', required: false)
curl = cc.find_library('curl')
# Cannot be found via pkg-config
ec = cc.find_library('ec')
ecm = cc.find_library('ecm')
# Cannot be found via pkg-config
glpk = cc.find_library('glpk')
ppl = cc.find_library('ppl')
gmpxx = cc.find_library('gmpxx')
readline = cc.find_library('readline')
fflas = dependency('fflas-ffpack')
fplll = dependency('fplll')
givaro = cc.find_library('givaro')
linbox = dependency('linbox', required: false)
if not linbox.found()
  linbox = cc.find_library('linbox')
endif
braiding = cc.find_library('braiding')
gc = cc.find_library('gc')
homfly = cc.find_library('homfly', has_headers: ['homfly.h'])
mpc = cc.find_library('mpc')
mpfi = cc.find_library('mpfi')
# Cannot be found via pkg-config (pkg-config file will be added in 4.13)
# Test for common.h header that was added in 4.12 as a indirect version check
gap = cc.find_library('gap', has_headers: ['gap/common.h'])

singular = dependency('Singular')

maxima = find_program('maxima', required: true)

# Cannot be found via pkg-config
ntl = cc.find_library('ntl', required: true)

# Meson currently ignores include_directories for Cython modules, so we
# have to add them manually.
# https://github.com/mesonbuild/meson/issues/9562
add_project_arguments('-I', meson.current_source_dir(), language: 'cython')
add_project_arguments('-I', meson.current_build_dir(), language: 'cython')

# Add global compiler flags
add_project_arguments('-X auto_pickle=False', language: 'cython')
add_project_arguments('-X autotestdict=False', language: 'cython')
add_project_arguments('-X binding=False', language: 'cython')
add_project_arguments('-X c_api_binop_methods=True', language: 'cython')
add_project_arguments('-X cdivision=True', language: 'cython')
add_project_arguments('-X cpow=True', language: 'cython')
add_project_arguments('-X embedsignature=True', language: 'cython')
add_project_arguments('--embed-positions', language: 'cython')
add_project_arguments('-X fast_getattr=True', language: 'cython')
#add_project_arguments('-X language_level="3"', language : 'cython')
add_project_arguments('-X legacy_implicit_noexcept=True', language: 'cython')
add_project_arguments(
  '-X preliminary_late_includes_cy28=True',
  language: 'cython',
)

inc_cpython = include_directories('sage/cpython')
inc_rings = include_directories('sage/rings')
inc_rings_finite = include_directories('sage/rings/finite_rings')
inc_flint = include_directories('sage/libs/flint')
inc_gsl = include_directories('sage/libs/gsl')
inc_ntl = include_directories('sage/libs/ntl')
inc_arb = include_directories('sage/libs/arb')
inc_data_structures = include_directories('sage/data_structures')
inc_ext = include_directories('sage/ext')
inc_partn_ref2 = include_directories('sage/groups/perm_gps/partn_ref2')
inc_src = include_directories('.')

# Submodules
subdir('sage')
