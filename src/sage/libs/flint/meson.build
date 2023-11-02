py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/flint',
)

extension_data = {
    'arith': files('arith.pyx'),
    'flint': files('flint.pyx'),
    'fmpz_factor': files('fmpz_factor.pyx'),
    'fmpz_poly': files('fmpz_poly.pyx'),
    'qsieve': files('qsieve.pyx'),
    'ulong_extras': files('ulong_extras.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/flint',
        install: true,
        include_directories: [inc_cpython, inc_src, inc_rings],
        dependencies: [py_dep, cysignals, gmp, flint],
    )
endforeach
