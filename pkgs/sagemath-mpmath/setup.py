#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-mpmath'],
           package_data={'sage.libs.mpmath': '*.pxd'})
