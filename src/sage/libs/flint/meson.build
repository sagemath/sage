py.install_sources(
    '__init__.py',
    'arith.pxd',
    'flint.pxd',
    'flint_ntl_wrap.h',
    'flint_wrap.h',
    'fmpq.pxd',
    'fmpq_mat.pxd',
    'fmpq_poly.pxd',
    'fmpz.pxd',
    'fmpz_factor.pxd',
    'fmpz_mat.pxd',
    'fmpz_mod.pxd',
    'fmpz_mod_poly.pxd',
    'fmpz_poly.pxd',
    'fmpz_poly_mat.pxd',
    'fmpz_poly_q.pxd',
    'fmpz_vec.pxd',
    'fq.pxd',
    'fq_nmod.pxd',
    'nmod_poly.pxd',
    'nmod_vec.pxd',
    'ntl_interface.pxd',
    'padic.pxd',
    'padic_poly.pxd',
    'qadic.pxd',
    'qsieve.pxd',
    'thread_pool.pxd',
    'types.pxd',
    'ulong_extras.pxd',
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
        include_directories: [inc_cpython, inc_flint, inc_rings],
        dependencies: [py_dep, cysignals, flint, gmp, mpfr],
    )
endforeach

