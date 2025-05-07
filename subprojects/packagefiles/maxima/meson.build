project('maxima', 'c',
  meson_version : '>=0.63.0',
)

mod = import('unstable-external_project')

cc = meson.get_compiler('c')
ecl = cc.find_library('ecl', required: false, disabler: true)
if not ecl.found()
  ecl_proj = subproject('ecl')
  ecl = ecl_proj.get_variable('ecl_dep')
endif
# It is not possible to have an external project directly depend on another external project (actually on an InternalDependency)
# So we need to create a dummy library as the middleman
ec_lib = library('ecl_lib',
  dependencies : [ecl.partial_dependency(sources : true)],
)
env = environment()
env.append('PATH', meson.current_build_dir() + '/../ecl/src/dist/usr/local/bin')
proj = mod.add_project('configure',
  configure_options : [
                        '--enable-ecl',
                        '--with-ecl=' + meson.current_build_dir() + '/../ecl/src/dist/usr/local/bin/ecl',
                        '--disable-build-docs',
                        # Maxima configure checks for git and, if it finds it, uses
                        # versions information from the repo. See #15529. We disable this.
                        'git_found=false',
                      ],
#  depends: [ec_lib],
#  env: env,
)

maxima_fas_dep = proj.dependency('maxima')
# It is not possible to have a custom_target depend on an external project (which returns a InternalDependency)
# So we need to create a dummy library as the middleman
#maxima_fas_lib = library('maxima_fas',
#  dependencies : [maxima_fas_dep.partial_dependency(sources : true)],
#)

ecl_bin = find_program('ecl', required : false)
if not ecl_bin.found()
  ecl_proj = subproject('ecl')
  ecl_bin = ecl_proj.get_variable('ecl_dep').get_variable('ecl_bin')
#else
#  # System ECL is not found, so we use a subproject to build it
#  ecldir = meson.current_build_dir() + '/../ecl/src/dist/usr/local/lib64/ecl-24.5.10'
endif
#ecldir = run_command([ecl_bin, '-eval', '(princ (SI:GET-LIBRARY-PATHNAME))', '-eval', '(quit)'], check: true).stdout().strip()

# Copy the Maxima fas file to the ECL directory
# TODO: Find a cleaner way to notify ecl about the maxima fas file
#copy_maxima_fas = custom_target('copy-maxima-fas',
#  output : 'maxima',
#  command : [
#    'cp',
#    meson.current_build_dir() + '/dist/usr/local/lib64/maxima/5.47.0/binary-ecl/maxima',
#    ecldir + '/maxima.fas',
#  ],
#  #input: [maxima_fas_lib],
#  depends: [maxima_fas_lib],
#  build_by_default : true,
#)
#copy_maxima_fas = custom_target('rename-maxima-fas',
#  output : 'maxima.fas',
#  command : [
#    'echo',
#    meson.current_build_dir() + '/dist/usr/local/lib64/maxima/5.47.0/binary-ecl/maxima',
#    meson.current_build_dir() + '/dist/usr/local/lib64/maxima/5.47.0/binary-ecl/maxima.fas',
#  ],
#  depends: [maxima_fas_lib],
#  build_by_default : true,
#)

#maxima_dep = declare_dependency(sources: [copy_maxima_fas])
#maxima_dep = declare_dependency(
#  variables: {
#    'MAXIMA_FAS' : meson.current_build_dir() + '/dist/usr/local/lib64/maxima/5.47.0/binary-ecl/maxima.fas',
#  },
#  sources: [copy_maxima_fas],
#)
maxima_dep = maxima_fas_dep.partial_dependency(sources : true)
meson.override_dependency('maxima', maxima_fas_dep)
