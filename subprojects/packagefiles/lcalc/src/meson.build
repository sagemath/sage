conf = configuration_data()
conf.set('PRECISION_DOUBLE', 1) # TODO: add more precision options
# TODO: add more data?

# Store this in "src" as opposed to "src/libLfunction" because
# doing the latter causes the other files in src/libLfunction
# to prefer the source-tree version of config.h from the dist
# tarball over the one in the build-tree that the user generates.
configure_file(
    output : 'config.h',
    configuration : conf,
    install_dir : get_option('includedir'),
    install : true)
config_h_dep = declare_dependency(
  include_directories : include_directories('.'))

subdir('libLfunction')
#subdir('lcalc')
