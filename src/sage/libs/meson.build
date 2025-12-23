sirocco = dependency('libsirocco', required: false, disabler: true)
if not sirocco.found()
  # fallback for older versions without pkg-config
  sirocco = cc.find_library(
    'sirocco',
    required: get_option('sirocco'),
    disabler: true,
  )
endif
# cannot be found via pkg-config
ecl = cc.find_library('ecl', required: false, disabler: true)
if not ecl.found() and not is_windows
  ecl_proj = subproject('ecl')
  ecl = ecl_proj.get_variable('ecl_dep')
endif
ecl_bin = find_program('ecl', required: false)
if ecl_bin.found()
  maxima_check = run_command(
    ecl_bin,
    '--eval',
    '(require \'maxima)',
    '--eval',
    '(quit)',
    check: false,
  )
  if maxima_check.returncode() == 0
    # ecl does not have any problems with Maxima, so nothing needs to be set here:
    conf_data.set('SAGE_MAXIMA_FAS', '')
  else
    # Try Sage's custom maxima location
    sage_maxima = join_paths(
      get_option('SAGE_LOCAL'),
      'lib',
      'ecl',
      'maxima.fas',
    )
    maxima_check = run_command(
      ecl_bin,
      '--eval',
      '(require \'maxima "' + sage_maxima + '")',
      '--eval',
      '(quit)',
      check: false,
    )
    conf_data.set('SAGE_MAXIMA_FAS', sage_maxima)
  endif

  maxima_bin = find_program('maxima', required: false, disabler: true)
  maxima_present = maxima_check.returncode() == 0 and maxima_bin.found()
else
  maxima_present = false
endif
if not maxima_present
  if is_windows
    # We currently cannot build Maxima on Windows, so we disable it
    maxima = disabler()
  else
    # ECL does not have Maxima support
    maxima_proj = subproject('maxima')
    maxima = maxima_proj.get_variable('maxima_dep')
    # Use the freshly built maxima
    conf_data.set('SAGE_MAXIMA_FAS', maxima.get_variable(internal: 'MAXIMA_FAS'))
    conf_data.set(
      'SAGE_MAXIMA_SHARE',
      maxima.get_variable(internal: 'MAXIMA_SHARE'),
    )
    conf_data.set('SAGE_MAXIMA', maxima.get_variable(internal: 'MAXIMA_BIN'))
  endif
else
  maxima = declare_dependency()
  # ecl does not have any problems with Maxima, so nothing needs to be set here:
  conf_data.set('SAGE_MAXIMA_SHARE', '')
  conf_data.set('SAGE_MAXIMA', maxima_bin.full_path())
endif

braiding = dependency(
  'libbraiding',
  version: '>=1.3.1',
  required: false,
  disabler: true,
)

# HiGHS optimization solver (used by sage.numerical.backends.highs_backend)
#
# Prefer pkg-config (highs.pc should provide -I$includedir/highs so relative
# includes like "lp_data/..." work). Some downstream distributions ship a
# broken highs.pc (e.g. pointing to /usr/debian/tmp), so we validate the
# pkg-config variables and fall back to finding the library directly.
highs = disabler()
highs_pc = dependency(
  'highs',
  required: false,
  method: 'pkg-config',
  disabler: true,
)
if highs_pc.found()
  highs_includedir = highs_pc.get_variable(
    pkgconfig: 'includedir',
    default_value: '',
  )
  highs_libdir = highs_pc.get_variable(pkgconfig: 'libdir', default_value: '')

  # Validate that the include directory exists and contains required headers
  highs_header_valid = (
    highs_includedir != '' and
    fs.is_dir(highs_includedir) and
    fs.is_file(highs_includedir / 'lp_data' / 'HighsCallbackStruct.h')
  )
  highs_lib_valid = highs_libdir != '' and fs.is_dir(highs_libdir)

  if highs_header_valid and highs_lib_valid
    highs = highs_pc
  else
    # pkg-config paths are invalid, fall back to finding the library directly
    highs_lib = cc.find_library('highs', required: false, disabler: true)
    if highs_lib.found()
      # Use the prefix from pkg-config or Sage's prefix to find headers
      highs_prefix = highs_pc.get_variable(pkgconfig: 'prefix', default_value: get_option('prefix'))
      highs_header_dir = highs_prefix / 'include' / 'highs'
      if fs.is_dir(highs_header_dir) and fs.is_file(highs_header_dir / 'lp_data' / 'HighsCallbackStruct.h')
        highs = declare_dependency(
          dependencies: [highs_lib],
          compile_args: ['-I' + highs_header_dir],
        )
      else
        highs = highs_lib
      endif
    endif
  endif
endif
gc = dependency(
  ['bdw-gc-threaded', 'bdw-gc'],
  version: '>=7.6.4',
  required: not is_windows,
  disabler: true,
)
homfly = dependency('libhomfly', required: false, disabler: true)
if not homfly.found()
  # fallback for older versions without pkg-config
  homfly = cc.find_library(
    'homfly',
    has_headers: ['homfly.h'],
    required: false,
    disabler: true,
  )
endif

py.install_sources(
  '__init__.py',
  'all.py',
  'braiding.pyx',
  'ecl.pxd',
  'ecl.pyx',
  'eclsig.h',
  'gmpxx.pxd',
  'homfly.pyx',
  'iml.pxd',
  'libecm.pyx',
  'm4ri.pxd',
  'm4rie.pxd',
  'meataxe.pxd',
  'meataxe.pyx',
  'sirocco.pyx',
  subdir: 'sage/libs',
)

extension_data = {
  'ecl' : files('ecl.pyx'),
  'homfly' : files('homfly.pyx'),
  'libecm' : files('libecm.pyx'),
  'meataxe': files('meataxe.pyx'),
}

dependencies = [py_dep, cysignals, gc, gmp]

foreach name, pyx : extension_data
  deps = dependencies
  if name == 'ecl'
    deps += [ecl, maxima]
  elif name == 'meataxe'
    deps += [mtx]
  elif name == 'homfly'
    deps += [homfly]
  elif name == 'libecm'
    deps += [ecm]
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: deps,
  )
endforeach

extension_data_cpp = {
  'braiding': files('braiding.pyx'),
  'sirocco': files('sirocco.pyx'),
}

foreach name, pyx : extension_data_cpp
  deps = dependencies
  if name == 'braiding'
    deps += [braiding]
  elif name == 'sirocco'
    deps += [sirocco]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: deps,
  )
endforeach

subdir('arb')
subdir('coxeter3')
install_subdir('cremona', install_dir: sage_install_dir / 'libs')
subdir('eclib')
subdir('flint')
subdir('gap')
install_subdir('giac', install_dir: sage_install_dir / 'libs')
subdir('glpk')
subdir('gmp')
subdir('gsl')
subdir('highs')
subdir('lcalc')
subdir('linbox')
install_subdir('linkages', install_dir: sage_install_dir / 'libs')
install_subdir('lrcalc', install_dir: sage_install_dir / 'libs')
install_subdir('mpfr', install_dir: sage_install_dir / 'libs')
subdir('mpmath')
install_subdir('mwrank', install_dir: sage_install_dir / 'libs')
subdir('ntl')
subdir('pari')
subdir('singular')
subdir('symmetrica')
