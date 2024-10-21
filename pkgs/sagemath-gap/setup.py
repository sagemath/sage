#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-gap'],
           package_data={'sage.libs.gap': ['sage.gaprc'],
                         'sage.ext_data.gap': ['*']})
