# -*- conf-unix -*-
[metadata]
name = sagemath-groups
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Groups and Invariant Theory
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

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
