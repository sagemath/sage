project('fflas-ffpack', 'c')

mod = import('unstable-external_project')

project = mod.add_project('configure')

fflasffpack_dep = project.dependency('fflas-ffpack')
meson.override_dependency('fflas-ffpack', fflasffpack_dep)
