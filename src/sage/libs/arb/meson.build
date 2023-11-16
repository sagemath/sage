py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/arb',
)

extension_data = {
    'arb_version': files('arb_version.pyx'),
    'arith': files('arith.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/arb',
        install: true,
        include_directories: [inc_arb, inc_cpython, inc_flint, inc_rings],
        dependencies: [py_dep, arb, flint, gmp, mpfr],
    )
endforeach

