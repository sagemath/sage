# Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cython = meson.get_compiler('cython')

# Setup dependencies that are needed by many modules
numpy = dependency('numpy', required: false, version: '>=2.2.4')
if not numpy.found()
  # Older versions of numpy don't have a pkg-config file
  inc_numpy = run_command(
    py,
    [
      '-c',
      '''
from os.path import relpath
import numpy
path = numpy.get_include()
try:
  print(relpath(path))
except Exception:
  print(path)
      '''.strip(),
    ],
    check: true,
  ).stdout().strip()
  numpy = declare_dependency(include_directories: inc_numpy)
endif

inc_cysignals = run_command(
  py,
  [
    '-c',
    '''
from os.path import relpath
import cysignals
path = cysignals.__file__.replace('__init__.py', '')
try:
  print(relpath(path))
except Exception:
  print(path)
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
from os.path import relpath
import gmpy2
path = gmpy2.__file__.replace('__init__.py', '')
try:
  print(relpath(path))
except Exception:
  print(path)
    '''.strip(),
  ],
  check: true,
).stdout().strip()
gmpy2 = declare_dependency(include_directories: inc_gmpy2)
gmp = dependency('gmp')

if is_windows
  # Not yet available on Windows
  cypari2 = disabler()
else
  inc_cypari2 = run_command(
    py,
    [
      '-c',
      '''
from os.path import relpath
import cypari2
path = cypari2.__file__.replace('__init__.py', '')
try:
  print(relpath(path))
except Exception:
  print(path)
      '''.strip(),
    ],
    check: true,
  ).stdout().strip()
  cypari2 = declare_dependency(include_directories: inc_cypari2)
endif
# Cannot be found via pkg-config
pari = cc.find_library('pari', required: not is_windows, disabler: true)
if pari.found() and not meson.is_cross_build()
  # Verify PARI version
  pari_version_code = '''
  #include <stdio.h>
  #include <pari/pari.h>
  int main(void) {
    pari_init(1000000, 2);
    GEN v = pari_version(), M = gel(v,1), m = gel(v,2), p = gel(v,3);
    printf("%ld.%ld.%ld", itos(M), itos(m), itos(p));
    pari_close();
    return 0;
  }
  '''
  pari_version = cc.run(
    pari_version_code,
    args: ['-v'],
    name: 'pari version',
    dependencies: [pari],
    required: true,
  ).stdout().strip()
  if pari_version.version_compare('<=2.17.0')
    message('PARI version > 2.17.0 required, found ' + pari_version)
    pari = disabler()
  endif
endif


mpfr = dependency('mpfr')

if is_windows
  # TODO: Reenable the following once conda's python is fixed
  # https://github.com/conda-forge/python-feedstock/pull/770
  # # In its currently released version, flint is not working on Windows; thus always use subproject
  # #flint = dependency('flint', version: '>=3.1.3')
  # cmake = import('cmake')
  # cmake_opts = cmake.subproject_options()
  # cmake_opts.add_cmake_defines({'BUILD_SHARED_LIBS': 'OFF'})
  # flint_proj = cmake.subproject('flint', options: cmake_opts)
  # flint = flint_proj.dependency('flint')
  # meson.override_dependency('flint', flint)
  flint = disabler()
else
  flint = dependency('flint', version: '>=3.0.0')
  if flint.version().version_compare('<3.1')
    # In older versions of flint, pkg-config file is broken, so we manually use find_library
    # This has been fixed in flint v3.1: https://github.com/flintlib/flint/pull/1647
    flint = cc.find_library('flint')
  endif
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
if is_windows
  # pkg-config file is wrong on Windows (https://github.com/conda-forge/gsl-feedstock/issues/63)
  gsl = cc.find_library('gsl', required: false, disabler: true)
else
  gsl = dependency('gsl', version: '>=2.5')
endif
gd = dependency('gdlib', required: false, version: '>=2.1')
if not gd.found()
  # Doesn't have a pkg-config file on some systems (https://github.com/conda-forge/libgd-feedstock/issues/55)
  gd = cc.find_library('gd', required: not is_windows, disabler: true)
endif
# Only some platforms have a standalone math library (https://mesonbuild.com/howtox.html#add-math-library-lm-portably)
m = cc.find_library('m', required: false)
m4ri = dependency('m4ri', version: '>=20250128', required: false)
if not m4ri.found()
  # Older versions of m4ri are declaring "@SIMD_CFLAGS@" as compile_args
  # which breaks the build. So we exclude these compile args.
  # https://github.com/malb/m4ri/commit/3197be5f7d3c67e8f13e32516557f35cbf4ff6ce
  m4ri = dependency('m4ri', version: '>=20140914').partial_dependency(
    includes: true,
    link_args: true,
    links: true,
    sources: true,
  )
endif
m4rie = dependency('m4rie', required: false)
if not m4rie.found()
  # For some reason, m4rie is not found via pkg-config on some systems (eg Conda)
  m4rie = cc.find_library('m4rie', required: not is_windows, disabler: true)
endif
mtx = dependency('libmtx', required: false)
if not mtx.found()
  # fallback for older versions without pkg-config
  mtx = cc.find_library(
    'mtx',
    required: get_option('meataxe'),
    disabler: true,
    has_headers: ['meataxe.h'],
  )
endif
png = dependency(['libpng', 'png', 'png16'], version: '>=1.2')
zlib = dependency('zlib', version: '>=1.2.11')
# We actually want >= 20250122, but the version number is not updated in the pkgconfig
# https://github.com/conda-forge/eclib-feedstock/issues/48
# To old versions of eclib don't ship a pkg-config file, so this should be rather safe
ec = dependency(
  'eclib',
  version: '>=20231211',
  required: get_option('eclib'),
  disabler: true,
)
ecm = cc.find_library('ecm', required: not is_windows, disabler: true)
gmpxx = dependency('gmpxx', required: not is_windows, disabler: true)
fflas = dependency(
  'fflas-ffpack',
  required: not is_windows,
  disabler: true,
  version: '>=2.5.0',
)
givaro = dependency(
  'givaro',
  required: not is_windows,
  disabler: true,
  version: '>=4.2.0',
)
linbox = dependency('linbox', required: false, version: '>=1.7.0')
if not linbox.found()
  linbox = cc.find_library('linbox', required: not is_windows, disabler: true)
endif
mpc = cc.find_library('mpc', required: not is_windows, disabler: true)
mpfi = cc.find_library('mpfi', required: false)
if not mpfi.found()
  mpfi_proj = subproject('mpfi')
  mpfi = mpfi_proj.get_variable('mpfi_dep')
endif

gap = dependency('libgap', version: '>=4.13.0', required: false)
if not gap.found()
  # Fallback in case pkg-config info is not available
  # Test for common.h header that was added in 4.12 as a indirect version check
  gap = cc.find_library(
    'gap',
    has_headers: ['gap/common.h'],
    required: not is_windows,
    disabler: true,
  )
endif
singular = dependency('Singular', required: not is_windows, disabler: true)
singular_factory = disabler()
if singular.found()
  singular_factory = singular
else
  singular_proj = subproject('singular')
  singular_factory = singular_proj.get_variable('factory_dep')
endif
# Cannot be found via pkg-config
ntl = cc.find_library('ntl', required: not is_windows, disabler: true)

boost_cpp17_compatible = false
boost = dependency(
  'boost',
  version: '>=1.83.0',
  required: false,
  modules: ['graph'],
)
if boost.found()
  # Need at least v1.83.0 (https://github.com/boostorg/functional/commit/6a573e4b8333ee63ee62ce95558c3667348db233)
  # for C++17 compatibility.
  boost_cpp17_compatible = true
else
  boost = dependency(
    'boost',
    version: '>=1.66.0',
    required: false,
    modules: ['graph'],
  )
endif


# Meson currently ignores include_directories for Cython modules, so we
# have to add them manually.
# https://github.com/mesonbuild/meson/issues/9562
add_project_arguments('-I', meson.current_source_dir(), language: 'cython')
add_project_arguments('-I', meson.current_build_dir(), language: 'cython')

# Add global compiler flags
add_project_arguments('-X auto_pickle=False', language: 'cython')
add_project_arguments('-X autotestdict=False', language: 'cython')
add_project_arguments('-X binding=True', language: 'cython')
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

src = meson.current_source_dir()

# Submodules
subdir('sage')
subdir('doc')
