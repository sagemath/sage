# -*- conf-unix -*-
[metadata]
name = sagemath-groups
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Groups
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
        sagemath_categories \
        sagemath_gap \
        sagemath_modules \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]
test = esyscmd(`sage-get-system-packages install-requires sagemath_repl')

# extras by packages
coxeter3 = esyscmd(`sage-get-system-packages install-requires sagemath_coxeter3')
gap =                            # no extra needed


# extras by groups_catalog
additive =                       # no extra needed
affine =                         # no extra needed
lie =  # FIXME
matrix =                         # no extra needed
permutation =                    # no extra needed
presentation =                   # no extra needed

# extras by other features
representations = esyscmd(`sage-get-system-packages install-requires sagemath_combinat')
semigroups = esyscmd(`sage-get-system-packages install-requires sagemath_combinat')

# the whole package
standard = sagemath-groups[additive,matrix,representations,semigroups]
