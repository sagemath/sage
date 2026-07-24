ecl_version = meson.project_version()
staged_prefix = meson.current_build_dir() + '/dist' + get_option('prefix')
libpath = staged_prefix / get_option('libdir')
ecldir = libpath / ('ecl-' + ecl_version)
# Keep side modules in the canonical SPKG/conda location, independently of the
# platform-specific libdir used by ECL itself.
ecl_side_modules_dir = get_option('prefix') / 'lib' / ('ecl-' + ecl_version)
include_path = meson.current_build_dir() / 'build'

ecl_env = environment()
ecl_env.append(
  'LDFLAGS',
  '-Wl,-rpath,' + libpath + ' -L' + libpath,
  separator : ' ',
)
# Workaround for https://gitlab.com/embeddable-common-lisp/ecl/-/issues/779
# Need to tell ECL where to find the fas files after installation.
ecl_env.append('CFLAGS', '-DECLDIR=\"' + ecldir + '\"', separator : ' ')
ecl_env.set('ECL_SIDE_MODULES_PATH', ecl_side_modules_dir)

# The external_project module only applies env to configure, not to the later
# make invocation where side-modules.lsp is compiled.  The wrapper persists
# ECL_SIDE_MODULES_PATH in the generated Makefile.
proj = mod.add_project('sage-ecl-configure',
  configure_options : [
                        # ecl embeds the include path, so we have to set it to the *installed* path
                        # workaround for https://gitlab.com/embeddable-common-lisp/ecl/-/issues/780
                        '--includedir=' + include_path,
                        '--enable-manual=no',
                        '--enable-unicode=yes',
                        '--with-defsystem',
                        '--with-extra-files=' + meson.current_source_dir() / 'util' / 'side-modules.lsp',
                       ],
  env : ecl_env,
)

ecl_dep = proj.dependency('ecl')
# Ugly workaround to include the rpath in the dependency
# https://github.com/mesonbuild/meson/issues/12970
ecl_rpath_dep = declare_dependency(
  link_args: '-Wl,-rpath,' + libpath + ' -L' + libpath,
  include_directories: [include_directories('build')],
)
ecl_dep = declare_dependency(
  dependencies : [ecl_dep, ecl_rpath_dep],
  variables: {
    'ecl_bin' : staged_prefix / get_option('bindir') / 'ecl',
    'ecldir' : ecldir,
    'ecl_version' : ecl_version,
  },
)
meson.override_dependency('ecl', ecl_dep)
