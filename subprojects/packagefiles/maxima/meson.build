project('maxima', 'c',
  meson_version : '>=1.7',
)

mod = import('unstable-external_project')

ecl_dep = dependency(
  'ecl',
  version: '>=26.5.5',
  allow_fallback: false,
)
ecl = ecl_dep.get_variable(internal: 'ecl_bin')
ecl_version = ecl_dep.version()

# An external project cannot directly depend on another external project's
# InternalDependency. Use a small Meson target to carry the build ordering.
ec_lib = library(
  'ecl_lib',
  'build_dependency.c',
  dependencies : [ecl_dep.partial_dependency(sources : true)],
)
dependencies = [ec_lib]

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
