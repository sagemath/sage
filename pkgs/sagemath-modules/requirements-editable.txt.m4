include(`sage_spkg_versions.m4')dnl
dnl Same as pyproject.toml.m4 build-system requires; FIXME: should pin to built wheels.
SPKG_INSTALL_REQUIRES_setuptools_wheel
SPKG_INSTALL_REQUIRES_wheel
SPKG_INSTALL_REQUIRES_pkgconfig
SPKG_INSTALL_REQUIRES_cython
SPKG_INSTALL_REQUIRES_gmpy2
SPKG_INSTALL_REQUIRES_cysignals
SPKG_INSTALL_REQUIRES_memory_allocator
-e ../sage-setup
dnl -e ../sagemath-environment
dnl -e ../sagemath-categories
dnl Same as setup.cfg.m4 install_requires
-e ../sage-conf
dnl SPKG_INSTALL_REQUIRES_gmpy2
SPKG_INSTALL_REQUIRES_cysignals
SPKG_INSTALL_REQUIRES_pplpy
dnl SPKG_INSTALL_REQUIRES_memory_allocator
dnl SPKG_INSTALL_REQUIRES_sagemath_environment
dnl -e ../sagemath-glpk    (gives AttributeError: 'Extension' object has no attribute '_needs_stub')
