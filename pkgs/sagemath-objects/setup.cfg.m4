# -*- conf-unix -*-
[metadata]
name = sagemath-objects
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: Sage objects, elements, parents, categories, coercion, metaclasses
long_description = file: README.rst
long_description_content_type = text/x-rst
include(`setup_cfg_metadata.m4')dnl'

[options]
include(`sage_spkg_versions.m4')dnl'
python_requires = >=3.9, <3.12
install_requires =
    SPKG_INSTALL_REQUIRES_gmpy2
    SPKG_INSTALL_REQUIRES_cysignals

[options.extras_require]
# Currently we do not use the sage doctester to test sagemath-objects,
# so we do not list sagemath-repl here.
test =


[options.package_data]
sage.cpython =
    pyx_visit.h
    string_impl.h
    cython_metaclass.h
    python_debug.h

sage.rings =
    integer_fake.h
