libpath = meson.current_build_dir() + '/dist/usr/local/lib64'
ecldir = meson.current_build_dir() + '/dist/usr/local/lib64/ecl-24.5.10'
include_path = meson.current_build_dir() + '/build'

proj = mod.add_project('configure',
  configure_options : [
                        # ecl embeds the include path, so we have to set it to the *installed* path
                        # workaround for https://gitlab.com/embeddable-common-lisp/ecl/-/issues/780
                        '--includedir=' + include_path,
                        '--enable-manual=no',
                        '--enable-unicode=yes',
                        '--with-defsystem',
                      ],
  env: {
    'LDFLAGS': '-Wl,-rpath,' + libpath + ' -L' + libpath,
    # Workaround for https://gitlab.com/embeddable-common-lisp/ecl/-/issues/779
    # Need to tell ECL where to find the fas files after installation
    'CFLAGS': '-DECLDIR=\"' + ecldir + '\"',
  }
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
    'ecl_bin' : meson.current_build_dir() + '/dist/usr/local/bin/ecl',
    'ecldir' : ecldir,
  },
)
meson.override_dependency('ecl', ecl_dep)
