sage_install_dir = py.get_install_dir() / 'sage'

# Generate the configuration file
conf_data = configuration_data()
conf_data.set('PACKAGE_VERSION', '1.2.3')
conf_data.set('SAGE_ROOT', meson.current_source_dir() / '..' / '..')
# We use Python's prefix here to make it work with conda
prefix = fs.as_posix(py.get_variable('prefix', ''))
conf_data.set('prefix', prefix)
datadir = fs.expanduser(get_option('datadir'))
if not fs.is_absolute(datadir)
  datadir = prefix / datadir
endif
conf_data.set('SAGE_SHARE', datadir)
if not fs.exists(datadir / 'cremona')
  message(
    'Warning: The specified datadir does not contain the necessary Cremona database. Either specify a different datadir or specify a correct the correct path via the environment variable SAGE_SHARE during runtime.',
  )
endif
conf_data.set('SAGE_MAXIMA', maxima.full_path())
# Conda's ecl does not have any problems with Maxima, so nothing needs to be set here:
conf_data.set('SAGE_MAXIMA_FAS', '')
# Kenzo cannot yet be provided by the system, so we always use the SAGE_LOCAL path for now.
conf_data.set('SAGE_KENZO_FAS', '\'${prefix}\'/lib/ecl/kenzo.fas')
# It can be found, so we don't have to set anything here:
conf_data.set('NTL_INCDIR', '')
conf_data.set('NTL_LIBDIR', '')
ecl_config = find_program('ecl-config', required: not is_windows, disabler: true)
conf_data.set('SAGE_ECL_CONFIG', ecl_config.full_path())
conf_data.set('SAGE_ARCHFLAGS', 'unset')
# not needed when using conda, as we then don't build any pc files
conf_data.set('SAGE_PKG_CONFIG_PATH', '')
openmp = dependency('openmp', required: false, disabler: true)
if openmp.found()
  conf_data.set('OPENMP_CFLAGS', '-fopenmp')
  conf_data.set('OPENMP_CXXFLAGS', '-fopenmp')
endif
gap_exe = find_program('gap', required: not is_windows, disabler: true)
if gap_exe.found()
  gaprun = run_command(
    gap_exe,
    '-r',
    '-q',
    '--bare',
    '--nointeract',
    '-c',
    'Display(JoinStringsWithSeparator(GAPInfo.RootPaths,";"));',
    check: true,
  )
  gap_root_paths = gaprun.stdout().strip()
  gap_root_paths = '${prefix}/lib/gap;${prefix}/share/gap;' + gaprun.stdout().strip()
  conf_data.set('GAP_ROOT_PATHS', gap_root_paths)
endif
ecm_bin = find_program(
  ['ecm', 'gmp-ecm'],
  required: not is_windows,
  disabler: true,
)
conf_data.set('SAGE_ECMBIN', ecm_bin.full_path())

config_file = configure_file(
  input: '../../pkgs/sage-conf/_sage_conf/_conf.py.in',
  output: 'config.py',
  install_dir: py.get_install_dir() / 'sage',
  install: true,
  configuration: conf_data,
)

# Packages that need no processing and can be installed directly
no_processing = [
  'databases',
  'doctest',
  'ext_data',
  'features',
  'game_theory',
  'homology',
  'knots',
  'logic',
  'manifolds',
  'parallel',
  'repl',
  'sandpiles',
  'tensor',
  'topology',
  'typeset',
]
foreach package : no_processing
  install_subdir(package, install_dir: sage_install_dir)
endforeach

py.install_sources(
  '__init__.py',
  'all.py',
  'all__sagemath_bliss.py',
  'all__sagemath_categories.py',
  'all__sagemath_coxeter3.py',
  'all__sagemath_environment.py',
  'all__sagemath_mcqd.py',
  'all__sagemath_meataxe.py',
  'all__sagemath_objects.py',
  'all__sagemath_repl.py',
  'all__sagemath_sirocco.py',
  'all__sagemath_tdlib.py',
  'all_cmdline.py',
  'env.py',
  'version.py',
  subdir: 'sage',
)

subdir('cpython')
subdir('libs')
subdir('misc')
subdir('structure')
subdir('algebras')
subdir('arith')
subdir('ext')
subdir('calculus')
subdir('categories')
subdir('coding')
subdir('combinat')
subdir('crypto')
subdir('data_structures')
subdir('functions')
subdir('games')
subdir('geometry')
subdir('graphs')
subdir('groups')
subdir('interacts')
subdir('interfaces')
subdir('lfunctions')
subdir('matrix')
subdir('matroids')
subdir('modular')
subdir('modules')
subdir('monoids')
subdir('numerical')
subdir('plot')
subdir('probability')
subdir('quadratic_forms')
subdir('quivers')
subdir('rings')
subdir('schemes')
subdir('sets')
subdir('stats')
subdir('symbolic')
subdir('tests')
subdir('dynamics')
subdir('sat')
