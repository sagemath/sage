#!/usr/bin/env python

from sage_setup import sage_setup

sage_setup(['sagemath-linbox'],
           required_modules=('fflas-ffpack', 'givaro', 'gsl', 'linbox', 'cblas',
                             'm4ri', 'gdlib', 'libpng', 'zlib'))
