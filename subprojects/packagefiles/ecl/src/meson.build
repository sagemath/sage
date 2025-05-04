proj = mod.add_project('configure',
  configure_options : [
                       '--enable-manual=no',
                       '--enable-unicode=yes',
                       '--with-defsystem' 
                      ],
)

ecl_dep = proj.dependency('ecl')
meson.override_dependency('ecl', ecl_dep)
