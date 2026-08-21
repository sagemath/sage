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

# Sage-the-distribution installs data below SAGE_LOCAL while Python lives in a
# nested virtual environment.  Check both locations before using the fallback.
database_datadirs = [datadir]
if get_option('SAGE_LOCAL') != ''
  database_datadirs += [fs.expanduser(get_option('SAGE_LOCAL')) / 'share']
endif
cremona_mini_found = false
ellcurves_found = false
graphs_found = false
polytopes_db_found = false
foreach database_datadir : database_datadirs
  cremona_mini_found = cremona_mini_found or fs.exists(
    database_datadir / 'cremona' / 'cremona_mini.db',
  )
  ellcurves_found = ellcurves_found or fs.exists(
    database_datadir / 'ellcurves' / 'rank0',
  )
  graphs_found = graphs_found or fs.exists(
    database_datadir / 'graphs' / 'graphs.db',
  )
  polytopes_db_found = polytopes_db_found or fs.exists(
    database_datadir / 'reflexive_polytopes' / 'Full3d',
  )
endforeach
# Keep fallback install targets in the main project: meson-python installs with
# --skip-subprojects (see pyproject.toml).  It maps 'datadir' to the Python
# data-scheme root; keep 'share' in the suffix, as for the Jupyter kernel below
# (meson-python issue #517).
meson_data_install_dir = get_option('datadir') / 'share'
if (
  not cremona_mini_found
  or not ellcurves_found
)
  elliptic_curves_project = subproject('elliptic_curves')
  elliptic_curves_generator = elliptic_curves_project.get_variable(
    'elliptic_curves_generator',
  )
  elliptic_curves_allcurves = elliptic_curves_project.get_variable(
    'elliptic_curves_allcurves',
  )
  custom_target(
    'cremona_mini_database',
    input: elliptic_curves_allcurves,
    output: 'cremona_mini.db',
    command: [
      py,
      elliptic_curves_generator,
      'cremona',
      '@INPUT@',
      '@OUTPUT@',
    ],
    install: true,
    install_dir: meson_data_install_dir / 'cremona',
  )
  custom_target(
    'elliptic_curves_rank_databases',
    input: [elliptic_curves_allcurves] + elliptic_curves_project.get_variable(
      'elliptic_curves_rank_sources',
    ),
    output: elliptic_curves_project.get_variable(
      'elliptic_curves_rank_outputs',
    ),
    command: [
      py,
      elliptic_curves_generator,
      'ranks',
      '@INPUT@',
      '--outputs',
      '@OUTPUT@',
    ],
    install: true,
    install_dir: meson_data_install_dir / 'ellcurves',
  )
endif
meson_graphs_data_dir = ''
if not graphs_found
  graphs_project = subproject('graphs')
  meson_graphs_data_dir = graphs_project.get_variable('graphs_data_dir')
  install_data(
    graphs_project.get_variable('graphs_files'),
    install_dir: meson_data_install_dir / 'graphs',
  )
endif
conf_data.set('MESON_GRAPHS_DATA_DIR', meson_graphs_data_dir)
meson_polytopes_db_data_dir = ''
if not polytopes_db_found
  polytopes_db_project = subproject('polytopes_db')
  meson_polytopes_db_data_dir = polytopes_db_project.get_variable(
    'polytopes_db_data_dir',
  )
  install_data(
    polytopes_db_project.get_variable('polytopes_db_files'),
    install_dir: meson_data_install_dir / 'reflexive_polytopes',
    preserve_path: true,
  )
endif
conf_data.set('MESON_POLYTOPES_DB_DATA_DIR', meson_polytopes_db_data_dir)
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
ecm_bin = find_program(
  ['ecm', 'gmp-ecm'],
  required: false,
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

# Record what build-time features (e.g. external package support) were
# enabled or disabled, and whether or not we should defer detection to
# runtime.
conf_data.set10('DEFER_FEATURE_CHECKS', get_option('defer_feature_checks'))
conf_data.set10('ECLIB_ENABLED', ec.found())

# Write config file
# Should be last so that subdir calls can modify the config data
config_file = configure_file(
  input: 'config.py.in',
  output: 'config.py',
  install_dir: sage_install_dir,
  install: true,
  configuration: conf_data,
)
