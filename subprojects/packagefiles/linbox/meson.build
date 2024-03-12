project('linbox', 'c')

mod = import('unstable-external_project')

project = mod.add_project('configure')

linbox_dep = project.dependency('linbox')
meson.override_dependency('linbox', linbox_dep)
