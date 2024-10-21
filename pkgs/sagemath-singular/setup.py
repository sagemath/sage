#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-singular'],
           required_modules=('Singular',
                             # from sagemath-linbox
                             'fflas-ffpack', 'givaro', 'gsl', 'linbox', 'cblas',
                             'm4ri', 'gdlib', 'libpng', 'zlib'))
