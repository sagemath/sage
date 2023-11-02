py.install_sources(
    '__init__.py',
    'function_factory.py',
    'standard_options.py',
    subdir: 'sage/libs/singular',
)

extension_data = {
    # 'function': files('function.pyx'), # not yet on conda
    # 'groebner_strategy': files('groebner_strategy.pyx'), # not yet on conda
    # 'option': files('option.pyx'), # not yet on conda
    # 'polynomial': files('polynomial.pyx'), # not yet on conda
    # 'ring': files('ring.pyx'), # not yet on conda
    # 'singular': files('singular.pyx'), # not yet on conda
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/singular',
        install: true,
        include_directories: [],
        dependencies: [py_dep, cysignals],
    )
endforeach

