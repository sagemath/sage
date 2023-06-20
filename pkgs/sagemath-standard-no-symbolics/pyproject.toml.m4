include(`sage_spkg_versions_toml.m4')dnl' -*- conf-toml -*-
[build-system]
# Minimum requirements for the build system to execute.
requires = [
    # Some version of sage-conf is required.
    # Note that PEP517/518 have no notion of optional sage_spkg dependencies:
    # https://github.com/pypa/pip/issues/6144
    SPKG_INSTALL_REQUIRES_sage_conf
    SPKG_INSTALL_REQUIRES_setuptools
    SPKG_INSTALL_REQUIRES_wheel
    SPKG_INSTALL_REQUIRES_sage_setup
    SPKG_INSTALL_REQUIRES_cypari
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_cython
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_jinja2
    SPKG_INSTALL_REQUIRES_jupyter_core
    SPKG_INSTALL_REQUIRES_numpy
    SPKG_INSTALL_REQUIRES_pkgconfig
    SPKG_INSTALL_REQUIRES_pplpy
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_environment
    SPKG_INSTALL_REQUIRES_sagemath_objects
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_modules
    SPKG_INSTALL_REQUIRES_sagemath_polyhedra
    SPKG_INSTALL_REQUIRES_sagemath_graphs
    SPKG_INSTALL_REQUIRES_sagemath_combinat
    SPKG_INSTALL_REQUIRES_sagemath_repl
    SPKG_INSTALL_REQUIRES_sagemath_schemes
]
build-backend = "setuptools.build_meta"
