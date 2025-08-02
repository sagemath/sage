project('maxima', 'c',
  meson_version : '>=0.63.0',
)

mod = import('unstable-external_project')

cc = meson.get_compiler('c')
ecl = find_program('ecl', required: false)
dependencies = []
if ecl.found()
  ecl = ecl.full_path()
else
  ecl_proj = subproject('ecl')
  ecl = ecl_proj.get_variable('ecl_dep').get_variable('ecl_bin')
  # It is not possible to have an external project directly depend on another external project (actually on an InternalDependency)
  # So we need to create a dummy library as the middleman
  ec_lib = library('ecl_lib',
    dependencies : [ecl_proj.get_variable('ecl_dep').partial_dependency(sources : true)],
  )
  dependencies += [ec_lib]
endif
proj = mod.add_project('configure',
  configure_options : [
                        '--enable-ecl',
                        '--with-ecl=' + ecl,
                        '--disable-build-docs',
                        # Maxima configure checks for git and, if it finds it, uses
                        # versions information from the repo. See #15529. We disable this.
                        'git_found=false',
                      ],
  depends: dependencies,
)

maxima_fas_dep = proj.dependency('maxima').partial_dependency(sources : true)
maxima_dep = declare_dependency(
  variables: {
    'MAXIMA_FAS': meson.current_build_dir() + '/build/src/binary-ecl/maxima.fas',
    'MAXIMA_SHARE': meson.current_build_dir() + '/dist/usr/local/share/maxima/5.47.0/share',
    'MAXIMA_BIN': meson.current_build_dir() + '/dist/usr/local/bin/maxima',
  },
  dependencies: maxima_fas_dep
)
meson.override_dependency('maxima', maxima_dep)
