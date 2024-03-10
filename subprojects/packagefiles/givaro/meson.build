project('givaro', 'c')

mod = import('unstable-external_project')

project = mod.add_project('configure')

givaro_dep = project.dependency('givaro')
meson.override_dependency('givaro', givaro_dep)
