#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-mpmath'],
           recurse_packages=['sage.libs.mpmath._vendor*'],
           package_data={'sage.libs.mpmath': ['*.pxd']})
