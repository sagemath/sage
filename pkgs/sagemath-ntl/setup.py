#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-ntl'],
           package_data={'sage.libs.ntl': ['*.h', '*.pxi']})
