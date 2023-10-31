py.install_sources(
    '__init__.py',
    'all.py',
    subdir: 'sage/libs/symmetrica',
)

extension_data = {
    'symmetrica': files('symmetrica.pyx'),
}

# Cannot be found by pkg-config
symmetrica = declare_dependency(
  dependencies: [
    cc.find_library('symmetrica'),
  ]
)

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/symmetrica',
        install: true,
        dependencies: [py_dep, cysignals, symmetrica],
    )
endforeach
