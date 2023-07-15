include(`sage_spkg_versions.m4')dnl' -*- conf-unix -*-
[metadata]
name = sagemath-modules
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Vectors, matrices, tensors, vector spaces, affine spaces, modules and algebras, additive groups, quadratic forms, homology, coding theory, matroids
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals
    SPKG_INSTALL_REQUIRES_memory_allocator
    SPKG_INSTALL_REQUIRES_sagemath_categories
    SPKG_INSTALL_REQUIRES_sagemath_mpmath

[options.extras_require]
test    = SPKG_INSTALL_REQUIRES_sagemath_repl

# extras by packages
flint   = SPKG_INSTALL_REQUIRES_sagemath_flint
gsl     = # No extra needed
linbox  = # FIXME: SPKG_INSTALL_REQUIRES_sagemath_linbox
m4ri    = # FIXME
m4rie   = # FIXME
meataxe = SPKG_INSTALL_REQUIRES_sagemath_meataxe
mpfi    = # FIXME
mpfr    = # No extra needed
mpmath  = # No extra needed
ntl     = # FIXME
numpy   = SPKG_INSTALL_REQUIRES_numpy
pari    = SPKG_INSTALL_REQUIRES_sagemath_pari

# extras by rings
RDF     = sagemath-modules[numpy]
CDF     = sagemath-modules[numpy]
RR      =                            # no extra needed
CC      =                            # no extra needed
RIF     =
CIF     =
RBF     = sagemath-modules[flint]
CBF     = sagemath-modules[flint]
GF      = sagemath-modules[pari]
GF2     = sagemath-modules[m4ri]
GF2e    = sagemath-modules[m4rie]
GF2n    = sagemath-modules[m4rie]
GFpn    = sagemath-modules[meataxe]
QQbar   = sagemath-modules[NumberField]
AA      = sagemath-modules[NumberField]
UCF     = sagemath-modules[NumberField]
Zp      = # FIXME
Qp      = sagemath-modules[Zp]
Zq      = sagemath-modules[Zp]
Qq      = sagemath-modules[Zp]
FiniteField     = sagemath-modules[GF]
NumberField     = # FIXME
QuadraticField  = sagemath-modules[NumberField]
CyclotomicField = sagemath-modules[NumberField]

# extras by features
invariant   = SPKG_INSTALL_REQUIRES_sagemath_groups
combinat    = SPKG_INSTALL_REQUIRES_sagemath_combinat
padics      = sagemath-modules[Zp]

# the whole package
standard    = sagemath-modules[invariant,combinat,padics,NumberField,FiniteField,m4ri,m4rie,flint,linbox,numpy,mpfi,ntl,pari]
