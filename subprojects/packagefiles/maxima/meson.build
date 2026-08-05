project('maxima', 'c',
  meson_version : '>=1.7',
)

mod = import('unstable-external_project')

cc = meson.get_compiler('c')
ecl = find_program('ecl', version: '>=26.5.5', required: false)
dependencies = []
if ecl.found()
  ecl_version = run_command(
    ecl,
    '--norc',
    '--eval',
    '(princ (lisp-implementation-version))',
    '--eval',
    '(quit)',
    check : true,
  ).stdout().strip()
  ecl = ecl.full_path()
else
  ecl_proj = subproject('ecl')
  ecl_dep = ecl_proj.get_variable('ecl_dep')
  ecl = ecl_dep.get_variable('ecl_bin')
  ecl_version = ecl_dep.get_variable('ecl_version')
  # It is not possible to have an external project directly depend on another external project (actually on an InternalDependency)
  # So we need to create a dummy library as the middleman
  ec_lib = library('ecl_lib', 'build_dependency.c',
    dependencies : [ecl_dep.partial_dependency(sources : true)],
  )
  dependencies += [ec_lib]
endif

staged_prefix = meson.current_build_dir() + '/dist' + get_option('prefix')
# Match the canonical SPKG/conda side-module location, not the platform libdir
# used by Maxima's regular image tree.
ecl_side_modules_subdir = 'lib' / ('ecl-' + ecl_version)
ecl_side_modules_dir = get_option('prefix') / ecl_side_modules_subdir

proj = mod.add_project('sage-maxima-configure',
  configure_options : [
                        '--enable-ecl',
                        '--with-ecl=' + ecl,
                        # Maxima configure checks for git and, if it finds it, uses
                        # versions information from the repo. See #15529. We disable this.
                        'git_found=false',
                      ],
  env : {
    'MAXIMA_FAS_DIR': ecl_side_modules_dir,
  },
  depends: dependencies,
)

maxima_fas_dep = proj.dependency('maxima').partial_dependency(sources : true)
maxima_dep = declare_dependency(
  variables: {
    'MAXIMA_FAS': staged_prefix / ecl_side_modules_subdir / 'maxima.fas',
    'MAXIMA_PREFIX': staged_prefix,
    'MAXIMA_BIN': staged_prefix / get_option('bindir') / 'maxima',
  },
  dependencies: maxima_fas_dep
)
meson.override_dependency('maxima', maxima_dep)

subdir('build')
