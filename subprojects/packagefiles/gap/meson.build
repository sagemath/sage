project('gap', 'c')

mod = import('unstable-external_project')

project = mod.add_project('configure')

gap_dep = project.dependency('gap')
meson.override_dependency('gap', gap_dep)
