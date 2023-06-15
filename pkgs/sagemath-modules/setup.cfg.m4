# -*- conf-unix -*-
[metadata]
name = sagemath-modules
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Vectors, matrices, tensors, vector spaces, affine spaces,
  modules and algebras, additive groups, quadratic forms, homology, coding theory, matroids
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
author = The Sage Developers
author_email = sage-support@googlegroups.com
url = https://www.sagemath.org

classifiers =
    Development Status :: 6 - Mature
    Intended Audience :: Education
    Intended Audience :: Science/Research
    License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: 3.11
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        gmpy2          \
        cysignals      \
        memory_allocator \
        sagemath_categories \
        mpmath \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

# extras by packages
flint = esyscmd(`sage-get-system-packages install-requires sagemath_flint')
linbox = esyscmd(`sage-get-system-packages install-requires sagemath_linbox')
m4ri = # FIXME
m4rie = # FIXME
meataxe = esyscmd(`sage-get-system-packages install-requires sagemath_meataxe')
mpfi = # FIXME
ntl = # FIXME
numpy = esyscmd(`sage-get-system-packages install-requires numpy')

# extras by rings
RDF = sagemath-modules[numpy]
CDF = sagemath-modules[numpy]
RR  =                            # no extra needed
CC  =                            # no extra needed
RIF =
CIF =
RBF = sagemath-modules[flint]
CBF = sagemath-modules[flint]
GF   = sagemath-modules[linbox]
GF2  = sagemath-modules[m4ri]
GF2e = sagemath-modules[m4rie]
GF2n = sagemath-modules[m4rie]
GFpn = sagemath-modules[meataxe]
FiniteField = sagemath-modules[GF]
NumberField = # FIXME
QuadraticField = sagemath-modules[NumberField]
QQbar = sagemath-modules[NumberField]
AA = sagemath-modules[NumberField]
CyclotomicField = sagemath-modules[NumberField]
UCF = sagemath-modules[NumberField]
Zp = # FIXME
Qp = sagemath-modules[Zp]
Zq = sagemath-modules[Zp]
Qq = sagemath-modules[Zp]

# extras by features
invariant = esyscmd(`sage-get-system-packages install-requires sagemath_groups')
combinat = esyscmd(`sage-get-system-packages install-requires sagemath_combinat')
padics = sagemath-modules[Zp]

# the whole package
standard = sagemath-modules[invariant,combinat,padics,NumberField,FiniteField,m4ri,m4rie,flint,linbox,numpy,mpfi,ntl]
