[build-system]
# Minimum requirements for the build system to execute.
requires = [
    "sage_setup[autogen]",
    # Some version of sage-conf is required.
    # Note that PEP517/518 have no notion of optional sage_spkg dependencies:
    # https://github.com/pypa/pip/issues/6144
     esyscmd(`sage-get-system-packages install-requires-toml \
        sage_conf      \
        setuptools \
        wheel          \
        sage_setup     \
        cypari         \
        cysignals      \
        cython         \
        gmpy2          \
        jinja2         \
        jupyter_core   \
        numpy          \
        pkgconfig      \
        pplpy          \
        memory_allocator \
                    ')]
build-backend = "setuptools.build_meta"

[tool.conda-lock]
platforms = [
    'osx-64', 'linux-64', 'linux-aarch64', 'osx-arm64'
]

[external]
# External dependencies in the format proposed by https://peps.python.org/pep-0725
# In the future, sage-the-distribution can read this information
build-requires = [
  "virtual:compiler/c",
  "virtual:compiler/cpp",
  "pkg:generic/pkg-config"
]

host-requires = [
  "virtual:interface/blas",
  "pkg:generic/boost",
  "pkg:generic/brial",
  "pkg:generic/cddlib",
  "pkg:generic/cliquer",
  "pkg:generic/ecl",
  "pkg:generic/eclib",
  "pkg:generic/ecm",
  "pkg:generic/fflas-ffpack",
  "pkg:generic/fplll",
  "pkg:generic/flint",
  "pkg:generic/libgd",
  "pkg:generic/gap",
  "pkg:generic/gfan",
  "pkg:generic/giac",
  "pkg:generic/givaro",
  "pkg:generic/glpk",
  "pkg:generic/gmp",
  "pkg:generic/gsl",
  "pkg:generic/iml",
  "pkg:generic/lcalc",
  "pkg:generic/libbraiding",
  "pkg:generic/libhomfly",
  "pkg:generic/linbox",
  "pkg:generic/lrcalc",
  "pkg:generic/m4ri",
  "pkg:generic/m4rie",
  "pkg:generic/maxima",
  "pkg:generic/mpc",
  "pkg:generic/mpfi",
  "pkg:generic/mpfr",
  "pkg:generic/nauty",
  "pkg:generic/ntl",
  "pkg:generic/palp",
  "pkg:generic/pari",
  "pkg:generic/pari-elldata",
  "pkg:generic/pari-galdata",
  "pkg:generic/pari-seadata",
  "pkg:generic/planarity",
  "pkg:generic/ppl",
  "pkg:generic/primesieve",
  "pkg:generic/primecount",
  "pkg:generic/qhull",
  "pkg:generic/rw",
  "pkg:generic/singular",
  "pkg:generic/symmetrica",
  "pkg:generic/sympow",
]

dependencies = [
  "pkg:generic/tachyon",
  "pkg:generic/sagemath-polytopes-db",
  "pkg:generic/sagemath-elliptic-curves",
]
