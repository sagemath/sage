proj = mod.add_project('configure',
  configure_options : [
                        '--enable-manual=no',
                        '--enable-unicode=yes',
                        '--with-defsystem',
                      ],
)

ecl_dep = proj.dependency('ecl')
# Ugly workaround to include the rpath in the dependency
# https://github.com/mesonbuild/meson/issues/12970
ecl_rpath_dep = declare_dependency(link_args: '-Wl,-rpath,' + meson.current_build_dir() + '/dist/usr/local/lib64')
ecl_dep = declare_dependency(
  dependencies : [ecl_dep, ecl_rpath_dep],
)
meson.override_dependency('ecl', ecl_dep)
