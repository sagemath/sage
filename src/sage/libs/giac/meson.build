py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/giac',
)


extension_data_cpp = {
    # 'giac': files('giac.pyx'), # not yet on conda
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/giac',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [],
        dependencies: [py_dep, cysignals],
    )
endforeach

