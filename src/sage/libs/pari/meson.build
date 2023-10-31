py.install_sources(
    '__init__.py',
    'all.py',
    'tests.py',
    subdir: 'sage/libs/pari',
)

extension_data = {
    'convert_flint': files('convert_flint.pyx'),
    'convert_gmp': files('convert_gmp.pyx'),
    'convert_sage': files('convert_sage.pyx'),
    'convert_sage_complex_double': files('convert_sage_complex_double.pyx'),
    'convert_sage_real_double': files('convert_sage_real_double.pyx'),
    'misc': files('misc.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/pari',
        install: true,
        include_directories: [inc_flint, inc_cpython],
        dependencies: [py_dep, cypari2, cysignals],
    )
endforeach
