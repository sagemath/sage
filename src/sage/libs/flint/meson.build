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
  'arith' : files('arith.pyx'),
  'arith_sage' : files('arith_sage.pyx'),
  'flint_sage' : files('flint_sage.pyx'),
  'fmpq_poly_sage' : files('fmpq_poly_sage.pyx'),
  'fmpz_factor_sage' : files('fmpz_factor_sage.pyx'),
  'fmpz_poly' : files('fmpz_poly.pyx'),
  'fmpz_poly_sage' : files('fmpz_poly_sage.pyx'),
  'qsieve' : files('qsieve.pyx'),
  'qsieve_sage' : files('qsieve_sage.pyx'),
  'ulong_extras' : files('ulong_extras.pyx'),
  'ulong_extras_sage' : files('ulong_extras_sage.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/flint',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings],
    dependencies: [py_dep, cysignals, flint, gmp, mpfr],
  )
endforeach

