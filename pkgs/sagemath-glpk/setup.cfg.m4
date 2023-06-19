# -*- conf-unix -*-
[metadata]
name = sagemath-glpk
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Linear and mixed integer linear optimization backend using GLPK
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
python_requires = >=3.8, <3.12
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        sagemath_objects \
        memory_allocator \
        cysignals \
        | sed "2,\$s/^/    /;"')dnl

[options.extras_require]

# No test requirements; see comment in tox.ini
test =

##     esyscmd(`sage-get-system-packages install-requires sagemath_repl')
##     esyscmd(`sage-get-system-packages install-requires sagemath_polyhedra')
