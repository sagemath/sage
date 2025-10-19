sage_install_dir = py.get_install_dir() / 'sage'

# Generate the configuration file
conf_data = configuration_data()
conf_data.set('EDITABLE_SRC', root)
conf_data.set('EDITABLE_BUILD', root_build)
conf_data.set('PACKAGE_VERSION', '1.2.3')
# We use Python's prefix here to make it work with conda
prefix = fs.as_posix(py.get_variable('prefix', ''))
conf_data.set('prefix', prefix)
datadir = fs.expanduser(get_option('datadir'))
if not fs.is_absolute(datadir)
  datadir = prefix / datadir
endif
if not fs.exists(datadir / 'cremona')
  message(
    'Warning: The specified datadir does not contain the necessary Cremona database. Either specify a different datadir or specify a correct the correct path via the environment variable SAGE_SHARE during runtime.',
  )
endif
# Kenzo cannot yet be provided by the system, so we always use the SAGE_LOCAL path for now.
conf_data.set('SAGE_KENZO_FAS', '\'${prefix}\'/lib/ecl/kenzo.fas')
# It can be found, so we don't have to set anything here:
conf_data.set('NTL_INCDIR', '')
conf_data.set('NTL_LIBDIR', '')
ecl_config = find_program('ecl-config', required: false)
if ecl_config.found()
  conf_data.set('SAGE_ECL_CONFIG', ecl_config.full_path())
endif
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
  conf_data.set('GAP_ROOT_PATHS', gap_root_paths)
endif
ecm_bin = find_program(
  ['ecm', 'gmp-ecm'],
  required: not is_windows,
  disabler: true,
)
conf_data.set('SAGE_ECMBIN', ecm_bin.full_path())
nauty_bins_prefix = ''
nauty_progs = [
  'directg',
  'gentourng',
  'geng',
  'genbg',
  'gentreeg',
  'genktreeg',
  'genposetg',
]
extra_bin_paths = []
if get_option('SAGE_LOCAL') != ''
  extra_bin_paths += [join_paths(get_option('SAGE_LOCAL'), 'bin')]
endif
foreach prog : nauty_progs
  nauty_check = find_program(prog, dirs: extra_bin_paths, required: false)
  if not nauty_check.found()
    nauty_check = find_program(
      'nauty-' + prog,
      dirs: extra_bin_paths,
      required: false,
    )
  endif
  if nauty_check.found()
    nauty_bins_prefix = nauty_check.full_path().replace(prog, '')
  endif
endforeach
conf_data.set('SAGE_NAUTY_BINS_PREFIX', nauty_bins_prefix)
# Packages that need no processing and can be installed directly
no_processing = [
  'cli',
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
  'all_cmdline.py',
  'all_test.py',
  'config_test.py',
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

# Jupyter kernel spec
# meson-python maps 'datadir' currently the python prefix, see
# https://github.com/mesonbuild/meson-python/issues/517
kernel_dir = get_option('datadir') / 'share/jupyter/kernels/sagemath'
install_data('ext_data/notebook-ipython/logo.svg', install_dir: kernel_dir)
install_data('ext_data/notebook-ipython/logo-64x64.png', install_dir: kernel_dir)
kernel_data = configuration_data()
kernel_data.set('SAGE_VERSION', meson.project_version())
configure_file(
  input: 'ext_data/notebook-ipython/kernel.json.in',
  output: 'kernel.json',
  install_dir: kernel_dir,
  install: true,
  configuration: kernel_data,
)

# Write config file
# Should be last so that subdir calls can modify the config data
config_file = configure_file(
  input: 'config.py.in',
  output: 'config.py',
  install_dir: sage_install_dir,
  install: true,
  configuration: conf_data,
)
