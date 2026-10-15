project(
  'eclib',
  'cpp',
  version: '20250627',
  license: 'GPL-2.0-or-later',
  default_options: ['cpp_std=c++17'],
)

cpp = meson.get_compiler('cpp')
flint = dependency('flint', version: '>=3.0.0')
if flint.version().version_compare('<3.1')
  # Match Sage's workaround for the older, broken pkg-config file.
  flint = cpp.find_library('flint')
endif
deps = [
  flint,
  dependency('gmp'),
  cpp.find_library('ntl'),
  cpp.find_library('pari'),
  dependency('threads'),
]

# The library sources from libsrc/Makefile.am. The other .cc files in that
# directory are templates included by these sources, not compilation units.
sources = files(
  'libsrc/interface.cc',
  'libsrc/unimod.cc',
  'libsrc/arith.cc',
  'libsrc/marith.cc',
  'libsrc/compproc.cc',
  'libsrc/vector.cc',
  'libsrc/matrix.cc',
  'libsrc/subspace.cc',
  'libsrc/kbessel.cc',
  'libsrc/svector.cc',
  'libsrc/smatrix.cc',
  'libsrc/smatrix_elim.cc',
  'libsrc/xsplit.cc',
  'libsrc/conic.cc',
  'libsrc/legendre.cc',
  'libsrc/quadratic.cc',
  'libsrc/illl.cc',
  'libsrc/hilbert.cc',
  'libsrc/timer.cc',
  'libsrc/cubic.cc',
  'libsrc/polys.cc',
  'libsrc/realroots.cc',
  'libsrc/parifact.cc',
  'libsrc/p2points.cc',
  'libsrc/gf.cc',
  'libsrc/xsplit_data.cc',
  'libsrc/threadpool.cc',
  'libsrc/logger.cc',
  'libsrc/curve.cc',
  'libsrc/curvedata.cc',
  'libsrc/curvered.cc',
  'libsrc/points.cc',
  'libsrc/cperiods.cc',
  'libsrc/isogs.cc',
  'libsrc/heights.cc',
  'libsrc/mwprocs.cc',
  'libsrc/lambda.cc',
  'libsrc/sifter.cc',
  'libsrc/sieve_search.cc',
  'libsrc/htconst.cc',
  'libsrc/egr.cc',
  'libsrc/saturate.cc',
  'libsrc/divpol.cc',
  'libsrc/pointsmod.cc',
  'libsrc/curvemod.cc',
  'libsrc/ffmod.cc',
  'libsrc/tlss.cc',
  'libsrc/elog.cc',
  'libsrc/getcurve.cc',
  'libsrc/egros.cc',
  'libsrc/mequiv.cc',
  'libsrc/mrank1.cc',
  'libsrc/mlocsol.cc',
  'libsrc/mglobsol.cc',
  'libsrc/mquartic.cc',
  'libsrc/mrank2.cc',
  'libsrc/qc.cc',
  'libsrc/sqfdiv.cc',
  'libsrc/version.cc',
  'libsrc/minim.cc',
  'libsrc/reduce.cc',
  'libsrc/transform.cc',
  'libsrc/desc2.cc',
  'libsrc/bitspace.cc',
  'libsrc/GetOpt.cc',
  'libsrc/twoadic.cc',
  'libsrc/descent.cc',
  'libsrc/newforms.cc',
  'libsrc/symb.cc',
  'libsrc/homspace.cc',
  'libsrc/cusp.cc',
  'libsrc/oldforms.cc',
  'libsrc/fixc6.cc',
  'libsrc/periods.cc',
  'libsrc/moddata.cc',
  'libsrc/pcprocs.cc',
  'libsrc/nfd.cc',
  'libsrc/curvesort.cc',
)

# As in the SPKG, use FLINT and NTL multiprecision, without Boost threading.
# Keep a shared library so all four Sage extensions share eclib's global state.
ec_lib = shared_library(
  'ec',
  sources,
  include_directories: include_directories('libsrc'),
  dependencies: deps,
  cpp_args: ['-DFLINT=1', '-DVERSION="' + meson.project_version() + '"'],
  install: true,
)

eclib_dep = declare_dependency(
  version: meson.project_version(),
  include_directories: include_directories('libsrc'),
  link_with: ec_lib,
  dependencies: deps,
)
meson.override_dependency('eclib', eclib_dep)

install_subdir('libsrc/eclib', install_dir: get_option('includedir'))
