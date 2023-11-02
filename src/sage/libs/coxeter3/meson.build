py.install_sources(
    '__init__.py',
    'all__sagemath_coxeter3.py',
    'coxeter_group.py',
    subdir: 'sage/libs/coxeter3',
)


extension_data_cpp = {
    # 'coxeter': files('coxeter.pyx'), # not yet on conda
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/coxeter3',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [],
        dependencies: [py_dep, cysignals, coxeter3],
    )
endforeach

