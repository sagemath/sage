include(`sage_spkg_versions.m4')dnl
dnl Same as setp.cfg.m4 install_requires
SPKG_INSTALL_REQUIRES_gmpy2
SPKG_INSTALL_REQUIRES_cysignals
SPKG_INSTALL_REQUIRES_memory_allocator
dnl To be added when ready for editable:
dnl -e ../sagemath-environment
dnl -e ../sagemath-objects
dnl -e ../sagemath-categories
dnl -e ../sagemath-mpmath
