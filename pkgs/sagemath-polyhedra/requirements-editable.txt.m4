include(`sage_spkg_versions.m4')dnl
dnl Same as setup.cfg.m4 install_requires; FIXME: should pin to built wheels.
SPKG_INSTALL_REQUIRES_gmpy2
SPKG_INSTALL_REQUIRES_cysignals
SPKG_INSTALL_REQUIRES_memory_allocator
SPKG_INSTALL_REQUIRES_pplpy
-e ../sagemath-modules
-e ../sage-conf
dnl -e ../sagemath-environment
dnl -e ../sagemath-objects
dnl -e ../sagemath-categories
dnl -e ../sagemath-mpmath
dnl -e ../sagemath-glpk    (gives AttributeError: 'Extension' object has no attribute '_needs_stub')
